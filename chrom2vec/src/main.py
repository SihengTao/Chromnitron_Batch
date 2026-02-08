# A python script that runs the ATAC-seq pipeline on replicates
import argparse
import os
import shutil
import subprocess
import time
import yaml

STAGE_DIRS = [
    "s1_FASTQ",
    "s2_fastp",
    "s3_hisat2",
    "s4_merged_bam",
    "s5_subsampled_bam",
    "s6_genrich",
    "s7_zarr",
    "s8_normalized",
    "s9_bigwig",
]

DELETE_STAGE_FLAGS = [
    ("delete_s1", "s1_FASTQ"),
    ("delete_s2", "s2_fastp"),
    ("delete_s3", "s3_hisat2"),
    ("delete_s4", "s4_merged_bam"),
    ("delete_s5", "s5_subsampled_bam"),
    ("delete_s6", "s6_genrich"),
    ("delete_s7", "s7_zarr"),
    ("delete_s8", "s8_normalized"),
    ("delete_s9", "s9_bigwig"),
]

def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Run the chrom2vec ATAC-seq pipeline.")
    parser.add_argument("config", help="Path to config YAML.")
    parser.add_argument("--skip-s9", action="store_true", help="Skip s9 bigwig generation.")
    parser.add_argument("--cleanup-after-s8", action="store_true",
                        help="Remove s1-s7 outputs after s8 completes.")
    parser.add_argument("--delete-fastq", action="store_true",
                        help="Delete input FASTQ files after s8 completes.")
    parser.add_argument("--delete-sra", action="store_true",
                        help="Delete input SRA files listed in cleanup_config.sra_files.")
    return parser.parse_args(argv)

def as_bool(value, default=False):
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in ("1", "true", "yes", "y")

def resolve_pipeline_options(config, skip_s9=None, cleanup_after_s8=None,
                             delete_fastq=None, delete_sra=None):
    pipeline_cfg = config.get("pipeline_config", {})
    cleanup_cfg = config.get("cleanup_config", {})

    keep_dirs = cleanup_cfg.get("keep_dirs")
    if keep_dirs is None:
        keep_dirs = ["s8_normalized"]

    use_flag_mode = any(flag in cleanup_cfg for flag, _ in DELETE_STAGE_FLAGS)
    if use_flag_mode:
        delete_dirs = [
            stage for flag, stage in DELETE_STAGE_FLAGS
            if as_bool(cleanup_cfg.get(flag), False)
        ]
        keep_dirs = [d for d in STAGE_DIRS if d not in set(delete_dirs)]
    else:
        delete_dirs = cleanup_cfg.get("delete_dirs")
        if delete_dirs is not None:
            delete_set = set(delete_dirs)
            delete_dirs = [d for d in STAGE_DIRS if d in delete_set]
            keep_dirs = [d for d in STAGE_DIRS if d not in delete_set]
        else:
            keep = set(keep_dirs)
            delete_dirs = [d for d in STAGE_DIRS if d not in keep]

    resolved = {
        "skip_s9": as_bool(pipeline_cfg.get("skip_s9", False)),
        "cleanup_after_s8": as_bool(cleanup_cfg.get("after_s8", False)),
        "delete_fastq": as_bool(cleanup_cfg.get("delete_fastq", False)),
        "delete_sra": as_bool(cleanup_cfg.get("delete_sra", False)),
        "keep_dirs": keep_dirs,
        "delete_dirs": delete_dirs,
        "copy_qc_to_s8": as_bool(cleanup_cfg.get("copy_qc_to_s8", False)),
        "qc_subdir": cleanup_cfg.get("qc_subdir") or "qc",
        "fastq_root": cleanup_cfg.get("fastq_root"),
        "sra_root": cleanup_cfg.get("sra_root"),
        "sra_files": cleanup_cfg.get("sra_files") or [],
    }

    if skip_s9 is not None:
        resolved["skip_s9"] = skip_s9
    if cleanup_after_s8 is not None:
        resolved["cleanup_after_s8"] = cleanup_after_s8
    if delete_fastq is not None:
        resolved["delete_fastq"] = delete_fastq
    if delete_sra is not None:
        resolved["delete_sra"] = delete_sra

    return resolved

def iter_fastq_paths(fastq_files):
    for rep_reads in (fastq_files or {}).values():
        if not isinstance(rep_reads, dict):
            continue
        for read_key in ("R1", "R2"):
            fastq_path = rep_reads.get(read_key)
            if fastq_path:
                yield fastq_path

def delete_paths(paths, allowed_root=None, label="file"):
    if allowed_root:
        allowed_root = os.path.realpath(allowed_root)
    for path in paths:
        if not path or not os.path.lexists(path):
            continue
        real_path = os.path.realpath(path)
        if allowed_root and os.path.commonpath([real_path, allowed_root]) != allowed_root:
            print(f"Skip deleting {label} outside allowed root: {path}")
            continue
        if os.path.isdir(path) and not os.path.islink(path):
            shutil.rmtree(path)
        else:
            os.remove(path)

def cleanup_intermediate_outputs(save_path, keep_dirs=None):
    keep = set(keep_dirs or [])
    for dir_name in STAGE_DIRS:
        if keep and dir_name in keep:
            continue
        target = os.path.join(save_path, dir_name)
        if os.path.exists(target):
            shutil.rmtree(target)
            print(f"Removed intermediate output: {target}")

def copy_qc_outputs_to_s8(save_path, qc_subdir="qc"):
    qc_root = os.path.join(save_path, "s8_normalized", qc_subdir)
    copied = 0

    fastp_root = os.path.join(save_path, "s2_fastp")
    if os.path.isdir(fastp_root):
        for rep in os.listdir(fastp_root):
            src_dir = os.path.join(fastp_root, rep)
            if not os.path.isdir(src_dir):
                continue
            dst_dir = os.path.join(qc_root, "s2_fastp", rep)
            os.makedirs(dst_dir, exist_ok=True)
            for name in ("fastp.json", "fastp.html"):
                src = os.path.join(src_dir, name)
                if os.path.isfile(src):
                    shutil.copy2(src, os.path.join(dst_dir, name))
                    copied += 1

    hisat2_root = os.path.join(save_path, "s3_hisat2")
    if os.path.isdir(hisat2_root):
        for rep in os.listdir(hisat2_root):
            src_dir = os.path.join(hisat2_root, rep)
            if not os.path.isdir(src_dir):
                continue
            src = os.path.join(src_dir, "summary.txt")
            if not os.path.isfile(src):
                continue
            dst_dir = os.path.join(qc_root, "s3_hisat2", rep)
            os.makedirs(dst_dir, exist_ok=True)
            shutil.copy2(src, os.path.join(dst_dir, "summary.txt"))
            copied += 1

    if copied:
        print(f"Copied {copied} QC files to {qc_root}")

def main():
    args = parse_args()
    config = load_config(args.config)
    run_pipeline(
        config,
        skip_s9=args.skip_s9,
        cleanup_after_s8=args.cleanup_after_s8,
        delete_fastq=args.delete_fastq,
        delete_sra=args.delete_sra,
    )

def run_pipeline(config, skip_s9=None, cleanup_after_s8=None,
                 delete_fastq=None, delete_sra=None):
    run_name = derive_run_name(config)
    # --- Begin execution with time stamp ---
    print(f"Starting pipeline at {time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}")

    fastq_files = config["pipeline_config"]["fastq_files"]
    save_path = f"{config['pipeline_config']['output_path']}/{run_name}"
    module_path = config["pipeline_config"]["module_path"]
    use_singularity = config["singularity_config"]["use_singularity"]
    bind_path = config["singularity_config"]["bind_path"]
    param_tuple = (fastq_files, save_path, module_path, use_singularity, bind_path, config)

    options = resolve_pipeline_options(
        config,
        skip_s9=skip_s9,
        cleanup_after_s8=cleanup_after_s8,
        delete_fastq=delete_fastq,
        delete_sra=delete_sra,
    )

    # Softlink, QC and alignment
    softlink_fastq_files(param_tuple)
    run_fastp(param_tuple)
    run_hisat2(param_tuple)

    # Merge and subsample bam files
    merge_bams(param_tuple)
    subsample_bam(param_tuple)

    # Run coverage and post-processing
    run_genrich(param_tuple)
    run_coverage_to_zarr(param_tuple)
    run_normalization(param_tuple)

    if not options["skip_s9"]:
        run_zarr_to_bigwig(param_tuple)

    if options["cleanup_after_s8"]:
        expected_s8 = os.path.join(
            save_path, "s8_normalized", "genrich_normalized.zarr"
        )
        if not os.path.exists(expected_s8):
            raise RuntimeError(f"s8 output not found; skipping cleanup: {expected_s8}")

        if options["copy_qc_to_s8"]:
            copy_qc_outputs_to_s8(save_path, options["qc_subdir"])

        keep_dirs = [d for d in STAGE_DIRS if d not in set(options["delete_dirs"])]
        cleanup_intermediate_outputs(save_path, keep_dirs=keep_dirs)

        if options["delete_fastq"]:
            fastq_paths = list(iter_fastq_paths(fastq_files))
            if not options["fastq_root"]:
                print("delete_fastq requested without fastq_root guard; deleting anyway.")
            delete_paths(
                fastq_paths,
                allowed_root=options["fastq_root"],
                label="FASTQ",
            )

        if options["delete_sra"]:
            if not options["sra_root"]:
                print("delete_sra requested without sra_root guard; deleting anyway.")
            delete_paths(
                options["sra_files"],
                allowed_root=options["sra_root"],
                label="SRA",
            )

# Decorator to calculate runtime for each function and terminal logging
def log_runtime(func):
    def wrapper(*args, **kwargs):
        print(f"\n{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}:\n--- Running {func.__name__} ---")
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        minutes = (end_time - start_time) / 60
        print(f"Function {func.__name__} took {minutes} minutes")
        return result
    return wrapper

def load_config(config_path):
    # Load config file
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def derive_run_name(config):
    """Return run_name; if set to 'auto' or empty, derive from FASTQ basenames."""
    run_name = str(config['pipeline_config'].get('run_name', '') or '').strip()
    if run_name and run_name.lower() != 'auto':
        config['pipeline_config']['run_name'] = run_name
        return run_name

    fastq_files = config['pipeline_config'].get('fastq_files') or {}
    candidates = []
    for rep_reads in fastq_files.values():
        for read_key in ('R1', 'R2'):
            fastq_path = rep_reads.get(read_key)
            if not fastq_path:
                continue
            base = os.path.basename(str(fastq_path))
            for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq', '.fastq.bz2', '.fq.bz2'):
                if base.endswith(ext):
                    base = base[:-len(ext)]
                    break
            for suffix in ('_1', '_R1', '-1', '-R1', '_001'):
                if base.endswith(suffix):
                    base = base[:-len(suffix)]
                    break
            candidates.append(base or 'run')
            break  # use first available read per replicate

    if not candidates:
        raise ValueError("run_name is empty/auto and no FASTQ files were found to derive it from.")

    first = candidates[0]
    if any(candidate != first for candidate in candidates):
        raise ValueError(f"Derived run_name differs across FASTQ files: {candidates}")

    config['pipeline_config']['run_name'] = first
    print(f"run_name not provided or set to auto; using derived value: {first}")
    return first

@log_runtime
def softlink_fastq_files(param_tuple):
    fastq_files, save_path, _, _, _, _ = param_tuple
    fastq_path = f'{save_path}/s1_FASTQ'
    # Link fastq files to new directory
    for rep, fastq_reads in fastq_files.items():
        os.makedirs(f"{fastq_path}/{rep}", exist_ok=True)
        for read, fastq_file in fastq_reads.items():
            os.makedirs(f"{fastq_path}/{rep}/{read}", exist_ok=True)
            dest = f"{fastq_path}/{rep}/{read}/fastq.gz"
            # If an old link/file exists, drop it unless it already points to the right source.
            if os.path.islink(dest):
                target = os.readlink(dest)
                target_abs = os.path.abspath(os.path.join(os.path.dirname(dest), target))
                src_abs = os.path.abspath(fastq_file)
                if target_abs == src_abs and os.path.exists(dest):
                    continue
                os.unlink(dest)
            elif os.path.exists(dest):
                os.remove(dest)
            if not os.path.exists(fastq_file):
                raise FileNotFoundError(f"FASTQ source not found: {fastq_file}")
            subprocess.run(["ln", "-s", f"{fastq_file}", dest]) # Soft link to fastq.gz file

@log_runtime
def run_fastp(param_tuple):
    fastq_files, save_path, module_path, use_singularity, bind_path, config = param_tuple
    fastq_path = f'{save_path}/s1_FASTQ'
    fastp_path = f'{save_path}/s2_fastp'

    for rep, fastq_reads in fastq_files.items():
        os.makedirs(f"{fastp_path}/{rep}", exist_ok=True)
        fastp_script_path = f"{module_path}/trimming/fastp_paired_end.sh"
        subprocess.run(["bash", fastp_script_path, use_singularity, bind_path,
                        f"{fastq_path}/{rep}/R1/fastq.gz", f"{fastq_path}/{rep}/R2/fastq.gz",
                        f"{fastp_path}/{rep}/R1.fastq.gz", f"{fastp_path}/{rep}/R2.fastq.gz",
                        f"{fastp_path}/{rep}/fastp.json", f"{fastp_path}/{rep}/fastp.html"])

@log_runtime
def run_hisat2(param_tuple):
    fastq_files, save_path, module_path, use_singularity, bind_path, config = param_tuple
    hisat2_path = f'{save_path}/s3_hisat2'
    fastp_path = f'{save_path}/s2_fastp'

    resources_path = config['pipeline_config']['resources_path']
    hisat2_index = config['resource_config']['hisat2_index']
    assembly_name = config['resource_config']['assembly_name']
    hisat2_index_path = f'{resources_path}/{hisat2_index}/{assembly_name}/genome'
    if not os.path.exists(f'{resources_path}/{hisat2_index}/{assembly_name}'):
        print('Index file does not exist, downloading...')
        download_hisat2_index(param_tuple)
    

    for rep, fastq_reads in fastq_files.items():
        os.makedirs(f"{hisat2_path}/{rep}", exist_ok=True)

        hisat2_script_path = f"{module_path}/alignment/hisat2_paired_end.sh"
        INDEX_PATH = hisat2_index_path
        subprocess.run(["bash", hisat2_script_path, use_singularity, bind_path, INDEX_PATH,
                        f"{fastp_path}/{rep}/R1.fastq.gz", f"{fastp_path}/{rep}/R2.fastq.gz",
                        f"{hisat2_path}/{rep}/hisat2.bam", f"{hisat2_path}/{rep}/summary.txt"])

@log_runtime
def download_hisat2_index(param_tuple):
    _, _, module_path, use_singularity, bind_path, config = param_tuple
    resources_path = config['pipeline_config']['resources_path']
    hisat2_index = config['resource_config']['hisat2_index']
    assembly_name = config['resource_config']['assembly_name']
    hisat2_index_root= f'{resources_path}/{hisat2_index}'
    hisat2_index_download_script_path = f"{module_path}/alignment/hisat2_download_index.sh"
    subprocess.run(["bash", hisat2_index_download_script_path, use_singularity, bind_path, hisat2_index_root, assembly_name])

@log_runtime
def merge_bams(param_tuple):
    fastq_files, save_path, module_path, use_singularity, bind_path, config = param_tuple
    hisat2_path = f'{save_path}/s3_hisat2'
    merged_bam_path = f'{save_path}/s4_merged_bam'
    os.makedirs(f"{merged_bam_path}", exist_ok=True)

    merge_script_path = f"{module_path}/bam/merge_bams.sh"
    inputs = [f"{hisat2_path}/{rep}/hisat2.bam" for rep in fastq_files.keys()]
    subprocess.run(["bash", merge_script_path, use_singularity, bind_path,
                    f"{merged_bam_path}/merged.bam", *inputs])

@log_runtime
def subsample_bam(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    merged_bam_path = f'{save_path}/s4_merged_bam'
    subsampled_bam_path = f'{save_path}/s5_subsampled_bam'
    os.makedirs(f"{subsampled_bam_path}", exist_ok=True)

    subsample_script_path = f"{module_path}/bam/subsample.sh"
    subsample_reads = "40000000"
    subprocess.run(["bash", subsample_script_path, use_singularity, bind_path,
                    "2023", subsample_reads,
                    f"{merged_bam_path}/merged.bam",
                    f"{subsampled_bam_path}/subsampled.bam"])

@log_runtime
def run_genrich(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    subsampled_bam_path = f'{save_path}/s5_subsampled_bam'
    genrich_path = f'{save_path}/s6_genrich'
    os.makedirs(f"{genrich_path}", exist_ok=True)

    genrich_script_path = f"{module_path}/coverage/genrich_coverage.sh"

    resources_path = config['pipeline_config']['resources_path']
    blacklist = config['resource_config']['blacklist']
    blacklist_path = f'{resources_path}/{blacklist}'
    assembly = config['resource_config']['assembly_name']

    os.makedirs(os.path.dirname(blacklist_path), exist_ok=True)
    subprocess.run(["bash", genrich_script_path, use_singularity, bind_path,
                    f"{subsampled_bam_path}/subsampled.bam", blacklist_path, assembly,
                    f"{genrich_path}/genrich.bedgraph"])

@log_runtime
def run_coverage_to_zarr(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    genrich_path = f'{save_path}/s6_genrich'
    zarr_path = f'{save_path}/s7_zarr'
    os.makedirs(f"{zarr_path}", exist_ok=True)

    coverage_script_path = f"{module_path}/io/coverage_to_zarr.sh"
    os.makedirs(f"{zarr_path}", exist_ok=True)
    subprocess.run(["bash", coverage_script_path, use_singularity, bind_path,
                    f"{genrich_path}/genrich.bedgraph",
                    f"{zarr_path}/genrich.zarr"])

@log_runtime
def run_normalization(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    zarr_path = f'{save_path}/s7_zarr'
    normalized_path = f'{save_path}/s8_normalized'
    os.makedirs(f"{normalized_path}", exist_ok=True)


    normalization_script_path = f"{module_path}/post_processing/normalization/normalize_atac_seq.sh"

    resources_path = config['pipeline_config']['resources_path']
    atac_peaks = config['resource_config']['atac_peaks']
    atac_non_peaks = config['resource_config']['atac_non_peaks']
    atac_peaks_path = f'{resources_path}/{atac_peaks}'
    atac_non_peaks_path = f'{resources_path}/{atac_non_peaks}'

    if not os.path.exists(atac_peaks_path):
        print('ATAC peak file does not exist, downloading...')
        peak_download_script_path = f"{module_path}/post_processing/normalization/download_atac_seq_peaks.sh"
        subprocess.run(["bash", peak_download_script_path, use_singularity, bind_path,
                    atac_peaks_path, atac_non_peaks_path])

    subprocess.run(["bash", normalization_script_path, use_singularity, bind_path,
                    atac_peaks_path, atac_non_peaks_path,
                    f"{zarr_path}/genrich.zarr",
                    f"{normalized_path}/genrich_normalized.zarr"])

@log_runtime
def run_zarr_to_bigwig(param_tuple):
    _, save_path, module_path, use_singularity, bind_path, config = param_tuple
    normalized_path = f'{save_path}/s8_normalized'
    bigwig_path = f'{save_path}/s9_bigwig'
    os.makedirs(f"{bigwig_path}", exist_ok=True)

    zarr_to_bigwig_script_path = f"{module_path}/io/zarr_to_bigwig.sh"  

    subprocess.run(["bash", zarr_to_bigwig_script_path, use_singularity, bind_path,
                    f"{normalized_path}/genrich_normalized.zarr",
                    f"{bigwig_path}/genrich_normalized.bw"])

if __name__ == "__main__":
    main()
