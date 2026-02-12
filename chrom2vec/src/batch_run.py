#!/usr/bin/env python3
import argparse
import copy
import glob
import os
import shlex
import subprocess
import sys
import time
import yaml

from concurrent.futures import ThreadPoolExecutor, as_completed

from main import run_pipeline, as_bool

DEFAULT_EBI_FTP_ROOT = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
DEFAULT_EBI_DOWNLOAD_CMD = ["wget"]
DEFAULT_EBI_DOWNLOAD_ARGS = ["-nc"]
DEFAULT_EBI_PROBE_CMD = ["wget"]
DEFAULT_EBI_PROBE_ARGS = ["--spider", "-q"]
DEFAULT_EBI_RETRY_COUNT = 2
DEFAULT_EBI_RETRY_WAIT_SEC = 30

DELETE_STAGE_FLAGS = [
    ("s1", "delete_s1"),
    ("s2", "delete_s2"),
    ("s3", "delete_s3"),
    ("s4", "delete_s4"),
    ("s5", "delete_s5"),
    ("s6", "delete_s6"),
    ("s7", "delete_s7"),
    ("s8", "delete_s8"),
    ("s9", "delete_s9"),
]

DELETE_STAGE_FLAG_KEYS = [flag for _, flag in DELETE_STAGE_FLAGS]

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Batch-run chrom2vec on SRR list.")
    parser.add_argument("config", help="Path to batch config YAML.")
    return parser.parse_args(argv)

def load_config(config_path):
    with open(config_path, "r") as f:
        return yaml.safe_load(f)

def get_cleanup_section(batch_cfg):
    cleanup = batch_cfg.get("cleanup")
    if isinstance(cleanup, dict):
        return cleanup
    return {}

def get_cleanup_value(batch_cfg, key, default=None):
    cleanup = get_cleanup_section(batch_cfg)
    if key in cleanup:
        return cleanup[key]
    if key == "after_s8" and "cleanup_after_s8" in batch_cfg:
        return batch_cfg["cleanup_after_s8"]
    if key in batch_cfg:
        return batch_cfg[key]
    return default

def get_stage_delete_value(batch_cfg, stage_short):
    cleanup = get_cleanup_section(batch_cfg)
    delete_stage = cleanup.get("delete_stage")
    if isinstance(delete_stage, dict) and stage_short in delete_stage:
        return delete_stage[stage_short]
    flat_key = f"delete_{stage_short}"
    if flat_key in batch_cfg:
        return batch_cfg[flat_key]
    return None

def should_use_stage_delete_flags(batch_cfg):
    cleanup = get_cleanup_section(batch_cfg)
    delete_stage = cleanup.get("delete_stage")
    if isinstance(delete_stage, dict):
        return True
    return any(key in batch_cfg for key in DELETE_STAGE_FLAG_KEYS)

def resolve_path(base_dir, path):
    if path is None:
        return None
    path = os.path.expanduser(os.path.expandvars(str(path)))
    if os.path.isabs(path):
        return path
    return os.path.abspath(os.path.join(base_dir, path))

def normalize_cmd(value, default_cmd):
    if value is None:
        return [default_cmd]
    if isinstance(value, list):
        return [str(item) for item in value]
    if isinstance(value, str):
        return shlex.split(value)
    raise ValueError(f"Unsupported command format: {value!r}")

def normalize_args(value):
    if not value:
        return []
    if isinstance(value, list):
        return [str(item) for item in value]
    if isinstance(value, str):
        return shlex.split(value)
    raise ValueError(f"Unsupported args format: {value!r}")

def build_download_cmd(batch_cfg, tool_key, default_cmd, bind_path=None):
    cmd = normalize_cmd(batch_cfg.get(tool_key), default_cmd)
    if as_bool(batch_cfg.get("download_use_singularity", False)):
        image = batch_cfg.get("download_singularity_image")
        if not image:
            raise ValueError("download_use_singularity requires download_singularity_image")
        singularity_cmd = normalize_cmd(batch_cfg.get("download_singularity_cmd"), "singularity")
        singularity_args = normalize_args(batch_cfg.get("download_singularity_args"))
        bind = bind_path or batch_cfg.get("download_singularity_bind")
        prefix = singularity_cmd + ["exec"] + singularity_args
        if bind:
            prefix += ["--bind", bind]
        cmd = prefix + [image] + cmd
    return cmd

def read_srr_list(path):
    srrs = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            srrs.append(line)
    return srrs

def find_sra_file(sra_root, srr):
    candidates = [
        os.path.join(sra_root, f"{srr}.sra"),
        os.path.join(sra_root, srr, f"{srr}.sra"),
    ]
    for path in candidates:
        if os.path.exists(path):
            return path

    matches = glob.glob(
        os.path.join(sra_root, "**", f"{srr}.sra"), recursive=True
    )
    if matches:
        return matches[0]

    raise FileNotFoundError(f"SRA file not found for {srr} under {sra_root}")

def run_command(cmd, label, capture_output=False):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    cmd_str = " ".join(cmd)
    print(f"{timestamp} | {label}: {cmd_str}")
    if capture_output:
        result = subprocess.run(cmd, text=True, capture_output=True)
        if result.stdout:
            print(result.stdout, end="")
        if result.stderr:
            print(result.stderr, end="")
        if result.returncode != 0:
            raise subprocess.CalledProcessError(
                result.returncode,
                cmd,
                output=result.stdout,
                stderr=result.stderr,
            )
        return result
    subprocess.run(cmd, check=True)

def is_unknown_arg_error(error_text, flag):
    if not error_text:
        return False
    if flag not in error_text:
        return False
    return "Unknown argument" in error_text or "rcParam" in error_text or "rcUnknown" in error_text

def compress_fastq_files(paths, compress_cmd, compress_args):
    for path in paths:
        if not os.path.exists(path):
            continue
        cmd = compress_cmd + compress_args + [path]
        run_command(cmd, f"compress {os.path.basename(path)}")
        if not os.path.exists(f"{path}.gz"):
            raise RuntimeError(f"Compression failed for {path}")

def run_prefetch(srr, sra_root, prefetch_cmd, prefetch_args, output_arg):
    cmd = prefetch_cmd + prefetch_args
    if output_arg:
        cmd += [output_arg, sra_root]
    cmd.append(srr)
    run_command(cmd, f"prefetch {srr}")

def run_fasterq_dump(srr, sra_path, output_dir, tmp_root, fasterq_cmd, fasterq_args,
                     threads, gzip_enabled, compress_cmd, compress_args):
    base_cmd = fasterq_cmd + fasterq_args + [
        "--split-files",
        "-e",
        str(threads),
        "-O",
        output_dir,
    ]
    if tmp_root:
        base_cmd += ["-t", tmp_root]

    cmd = list(base_cmd)
    if gzip_enabled:
        cmd.append("--gzip")
    cmd.append(sra_path)

    if gzip_enabled:
        try:
            run_command(cmd, f"fasterq-dump {srr}", capture_output=True)
            return True
        except subprocess.CalledProcessError as exc:
            error_text = (exc.stderr or "") + (exc.output or "")
            if is_unknown_arg_error(error_text, "--gzip"):
                print("fasterq-dump does not support --gzip; retrying without gzip.")
            else:
                raise

    cmd = list(base_cmd) + [sra_path]
    run_command(cmd, f"fasterq-dump {srr} (no gzip)")
    return False

def expected_fastq_paths(fastq_root, srr):
    out_dir = os.path.join(fastq_root, srr)
    r1 = os.path.join(out_dir, f"{srr}_1.fastq.gz")
    r2 = os.path.join(out_dir, f"{srr}_2.fastq.gz")
    return out_dir, r1, r2

def ebi_base_url(srr, ebi_ftp_root=DEFAULT_EBI_FTP_ROOT):
    if not str(srr).startswith("SRR"):
        return None
    num = str(srr)[3:]
    if not num.isdigit():
        return None

    prefix = num[:3]
    if len(num) <= 6:
        suffix = num[3:6]
    else:
        suffix = num[6:]
    if not suffix:
        return None

    try:
        suffix = f"{int(suffix):03d}"
    except ValueError:
        return None

    root = str(ebi_ftp_root).rstrip("/")
    return f"{root}/SRR{prefix}/{suffix}/{srr}"

def ebi_fastq_urls(srr, ebi_ftp_root=DEFAULT_EBI_FTP_ROOT):
    base_url = ebi_base_url(srr, ebi_ftp_root=ebi_ftp_root)
    if not base_url:
        return None, None
    return (
        f"{base_url}/{srr}_1.fastq.gz",
        f"{base_url}/{srr}_2.fastq.gz",
    )

def probe_remote_file(url, probe_cmd=None, probe_args=None):
    probe_cmd = probe_cmd or DEFAULT_EBI_PROBE_CMD
    probe_args = probe_args or DEFAULT_EBI_PROBE_ARGS
    cmd = probe_cmd + probe_args + [url]
    result = subprocess.run(cmd, text=True, capture_output=True)
    return result.returncode == 0

def try_download_from_ebi(srr, fastq_root):
    out_dir, r1, r2 = expected_fastq_paths(fastq_root, srr)
    if os.path.exists(r1) and os.path.exists(r2):
        return True

    r1_url, r2_url = ebi_fastq_urls(srr)
    if not (r1_url and r2_url):
        print(f"Skipping EBI check for invalid SRR format: {srr}")
        return False

    print(f"Checking EBI FASTQ availability for {srr}")
    if not (probe_remote_file(r1_url) and probe_remote_file(r2_url)):
        print(f"EBI FASTQ not found for {srr}, will retry or fallback.")
        return False

    os.makedirs(out_dir, exist_ok=True)
    for url in (r1_url, r2_url):
        cmd = DEFAULT_EBI_DOWNLOAD_CMD + DEFAULT_EBI_DOWNLOAD_ARGS + ["-P", out_dir, url]
        run_command(cmd, f"ebi download {srr}")

    return os.path.exists(r1) and os.path.exists(r2)

def download_srr(srr, sra_root, fastq_root, tmp_root,
                 prefetch_cmd, prefetch_args, prefetch_output_arg,
                 fasterq_cmd, fasterq_args, threads, gzip_enabled,
                 compress_cmd, compress_args):
    os.makedirs(sra_root, exist_ok=True)
    os.makedirs(fastq_root, exist_ok=True)

    out_dir, r1, r2 = expected_fastq_paths(fastq_root, srr)
    if os.path.exists(r1) and os.path.exists(r2):
        return r1, r2, None

    for attempt in range(DEFAULT_EBI_RETRY_COUNT + 1):
        try:
            if try_download_from_ebi(
                srr=srr,
                fastq_root=fastq_root,
            ):
                print(f"Using EBI FASTQ for {srr}")
                return r1, r2, None
        except Exception as exc:
            print(f"EBI attempt {attempt + 1} failed for {srr}: {exc}")

        if attempt < DEFAULT_EBI_RETRY_COUNT:
            print(
                f"Retrying EBI for {srr} in {DEFAULT_EBI_RETRY_WAIT_SEC}s "
                f"({attempt + 1}/{DEFAULT_EBI_RETRY_COUNT})"
            )
            time.sleep(DEFAULT_EBI_RETRY_WAIT_SEC)
        else:
            print(f"EBI failed after retries for {srr}, fallback to NCBI.")

    # Avoid partially downloaded pairs interfering with NCBI fallback.
    if os.path.exists(r1) != os.path.exists(r2):
        for path in (r1, r2):
            if os.path.exists(path):
                os.remove(path)

    sra_path = None
    try:
        sra_path = find_sra_file(sra_root, srr)
    except FileNotFoundError:
        pass

    if not sra_path:
        run_prefetch(
            srr,
            sra_root=sra_root,
            prefetch_cmd=prefetch_cmd,
            prefetch_args=prefetch_args,
            output_arg=prefetch_output_arg,
        )
        sra_path = find_sra_file(sra_root, srr)

    if not (os.path.exists(r1) and os.path.exists(r2)):
        os.makedirs(out_dir, exist_ok=True)
        used_gzip = run_fasterq_dump(
            srr=srr,
            sra_path=sra_path,
            output_dir=out_dir,
            tmp_root=tmp_root,
            fasterq_cmd=fasterq_cmd,
            fasterq_args=fasterq_args,
            threads=threads,
            gzip_enabled=gzip_enabled,
            compress_cmd=compress_cmd,
            compress_args=compress_args,
        )
        if not used_gzip:
            r1_raw = r1[:-3]
            r2_raw = r2[:-3]
            compress_fastq_files([r1_raw, r2_raw], compress_cmd, compress_args)

    if not (os.path.exists(r1) and os.path.exists(r2)):
        raise FileNotFoundError(f"FASTQ outputs missing for {srr}: {r1}, {r2}")

    return r1, r2, sra_path

def build_sample_config(base_config, srr, r1, r2, batch_cfg, sra_root, fastq_root, sra_path):
    sample_config = copy.deepcopy(base_config)
    sample_config.setdefault("pipeline_config", {})
    sample_config["pipeline_config"]["run_name"] = srr
    sample_config["pipeline_config"]["fastq_files"] = {
        "rep1": {"R1": r1, "R2": r2}
    }

    run_s9 = batch_cfg.get("run_s9")
    if run_s9 is not None:
        sample_config["pipeline_config"]["skip_s9"] = not as_bool(run_s9, False)

    cleanup_cfg = sample_config.setdefault("cleanup_config", {})
    cleanup_cfg["after_s8"] = as_bool(
        get_cleanup_value(batch_cfg, "after_s8", True)
    )

    use_flag_mode = should_use_stage_delete_flags(batch_cfg)
    if use_flag_mode:
        for stage_short, flag_key in DELETE_STAGE_FLAGS:
            value = get_stage_delete_value(batch_cfg, stage_short)
            cleanup_cfg[flag_key] = as_bool(value, False)
        cleanup_cfg.pop("keep_dirs", None)
        cleanup_cfg.pop("delete_dirs", None)
    else:
        delete_dirs = get_cleanup_value(batch_cfg, "delete_dirs")
        if delete_dirs is not None:
            cleanup_cfg["delete_dirs"] = delete_dirs
            cleanup_cfg.pop("keep_dirs", None)
        else:
            cleanup_cfg.pop("delete_dirs", None)
            cleanup_cfg["keep_dirs"] = (
                get_cleanup_value(batch_cfg, "keep_dirs") or ["s8_normalized"]
            )

    cleanup_cfg["delete_fastq"] = as_bool(
        get_cleanup_value(batch_cfg, "delete_fastq", True)
    )
    cleanup_cfg["delete_sra"] = as_bool(
        get_cleanup_value(batch_cfg, "delete_sra", True)
    )
    cleanup_cfg["copy_qc_to_s8"] = as_bool(
        get_cleanup_value(batch_cfg, "copy_qc_to_s8", False)
    )
    qc_subdir = get_cleanup_value(batch_cfg, "qc_subdir")
    if qc_subdir:
        cleanup_cfg["qc_subdir"] = str(qc_subdir)
    cleanup_cfg["fastq_root"] = fastq_root
    cleanup_cfg["sra_root"] = sra_root
    cleanup_cfg["sra_files"] = [sra_path] if sra_path else []
    return sample_config

def should_run_s9(base_config, batch_cfg):
    run_s9 = batch_cfg.get("run_s9")
    if run_s9 is not None:
        return as_bool(run_s9, False)
    skip_s9 = base_config.get("pipeline_config", {}).get("skip_s9")
    return not as_bool(skip_s9, False)

def should_keep_stage_after_cleanup(batch_cfg, stage_short, stage_dir):
    use_flag_mode = should_use_stage_delete_flags(batch_cfg)
    if use_flag_mode:
        return not as_bool(get_stage_delete_value(batch_cfg, stage_short), False)

    delete_dirs = get_cleanup_value(batch_cfg, "delete_dirs")
    if delete_dirs is not None:
        return stage_dir not in set(delete_dirs)

    keep_dirs = get_cleanup_value(batch_cfg, "keep_dirs") or ["s8_normalized"]
    return stage_dir in set(keep_dirs)

def get_completion_marker_path(output_root, srr, base_config, batch_cfg):
    run_s9 = should_run_s9(base_config, batch_cfg)
    cleanup_after_s8 = as_bool(get_cleanup_value(batch_cfg, "after_s8", True))

    keep_s8 = True
    keep_s9 = run_s9
    if cleanup_after_s8:
        keep_s8 = should_keep_stage_after_cleanup(batch_cfg, "s8", "s8_normalized")
        if run_s9:
            keep_s9 = should_keep_stage_after_cleanup(batch_cfg, "s9", "s9_bigwig")

    if run_s9 and keep_s9:
        return os.path.join(
            output_root,
            srr,
            "s9_bigwig",
            "genrich_normalized.bw",
        )
    if keep_s8:
        return os.path.join(
            output_root,
            srr,
            "s8_normalized",
            "genrich_normalized.zarr",
        )
    return None

def process_sample(srr, base_config, batch_cfg, config_dir):
    output_root = base_config.get("pipeline_config", {}).get("output_path")
    if output_root:
        expected = get_completion_marker_path(output_root, srr, base_config, batch_cfg)
        if expected and os.path.exists(expected):
            print(f"Skipping {srr}: existing output found at {expected}")
            return srr

    sra_root = resolve_path(config_dir, batch_cfg.get("sra_root", "sra"))
    fastq_root = resolve_path(config_dir, batch_cfg.get("fastq_root", "fastq"))
    tmp_root = resolve_path(config_dir, batch_cfg.get("temp_root"))
    if tmp_root:
        tmp_root = os.path.join(tmp_root, srr)
        os.makedirs(tmp_root, exist_ok=True)
    download_bind = batch_cfg.get("download_singularity_bind")
    if not download_bind:
        download_bind = base_config.get("singularity_config", {}).get("bind_path")
    prefetch_cmd = build_download_cmd(batch_cfg, "prefetch_cmd", "prefetch", bind_path=download_bind)
    fasterq_cmd = build_download_cmd(batch_cfg, "fasterq_dump_cmd", "fasterq-dump", bind_path=download_bind)
    prefetch_output_arg = batch_cfg.get("prefetch_output_arg", "--output-directory")
    prefetch_args = normalize_args(batch_cfg.get("prefetch_args"))
    fasterq_args = [arg for arg in normalize_args(batch_cfg.get("fasterq_dump_args")) if arg != "--gzip"]
    threads = int(batch_cfg.get("fasterq_threads", 8))
    gzip_enabled = as_bool(batch_cfg.get("fasterq_use_gzip", True))
    compress_cmd = normalize_cmd(batch_cfg.get("compress_cmd"), "gzip")
    compress_args = normalize_args(batch_cfg.get("compress_args"))

    r1, r2, sra_path = download_srr(
        srr,
        sra_root=sra_root,
        fastq_root=fastq_root,
        tmp_root=tmp_root,
        prefetch_cmd=prefetch_cmd,
        prefetch_args=prefetch_args,
        prefetch_output_arg=prefetch_output_arg,
        fasterq_cmd=fasterq_cmd,
        fasterq_args=fasterq_args,
        threads=threads,
        gzip_enabled=gzip_enabled,
        compress_cmd=compress_cmd,
        compress_args=compress_args,
    )

    sample_config = build_sample_config(
        base_config,
        srr=srr,
        r1=r1,
        r2=r2,
        batch_cfg=batch_cfg,
        sra_root=sra_root,
        fastq_root=fastq_root,
        sra_path=sra_path,
    )

    run_pipeline(sample_config)

    return srr

def main():
    args = parse_args()
    config_path = os.path.abspath(args.config)
    config_dir = os.path.dirname(config_path)
    config = load_config(config_path)

    batch_cfg = config.get("batch_config", {})
    srr_list_path = resolve_path(
        config_dir, batch_cfg.get("srr_list_path", "srr_list.txt")
    )
    srrs = read_srr_list(srr_list_path)
    if not srrs:
        raise ValueError(f"No SRR entries found in {srr_list_path}")
    if len(set(srrs)) != len(srrs):
        raise ValueError("Duplicate SRR IDs found in srr_list.txt")

    base_config = {k: v for k, v in config.items() if k != "batch_config"}
    pipeline_cfg = base_config.get("pipeline_config", {})
    for key in ("output_path", "module_path", "resources_path"):
        if key in pipeline_cfg:
            pipeline_cfg[key] = resolve_path(config_dir, pipeline_cfg[key])

    workers = int(batch_cfg.get("parallel_workers", 1))
    workers = max(1, min(workers, len(srrs)))

    print(f"Starting batch run with {workers} workers for {len(srrs)} samples.")

    results = []
    failures = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        future_map = {
            executor.submit(process_sample, srr, base_config, batch_cfg, config_dir): srr
            for srr in srrs
        }
        for future in as_completed(future_map):
            srr = future_map[future]
            try:
                result = future.result()
                results.append(result)
                print(f"Completed {result}")
            except Exception as exc:
                failures.append((srr, str(exc)))
                print(f"Failed {srr}: {exc}")

    if failures:
        print("Failures:")
        for srr, err in failures:
            print(f"  {srr}: {err}")
        sys.exit(1)

if __name__ == "__main__":
    main()
