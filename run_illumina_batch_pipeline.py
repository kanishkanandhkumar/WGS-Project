#!/usr/bin/env python3
import os
import sys
import subprocess
import logging
import shutil
import time
import argparse
import glob

# ==========================================
# CONSTANTS & DIRECTORY CONFIGURATION
# ==========================================
DIRS = ["reference", "alignment", "variants", "results", "reports", "logs", "results/batch_outputs"]

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler("logs/illumina_vep_pipeline.log"),
            logging.StreamHandler(sys.stdout)
        ]
    )

def run_cmd(cmd, step_name):
    start_time = time.time()
    logging.info(f"STARTING: {step_name}")
    try:
        process = subprocess.run(
            f"set -o pipefail; {cmd}", 
            shell=True, check=True, executable='/bin/bash',
            stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        duration = time.time() - start_time
        logging.info(f"SUCCESS: {step_name} completed in {duration:.2f}s")
    except subprocess.CalledProcessError as e:
        logging.error(f"ERROR: Step failed during {step_name}")
        logging.error(f"Exit code: {e.returncode}")
        logging.error(f"Error output:\n{e.stderr.decode('utf-8')}")
        sys.exit(1)

def prepare_vep_references(fna_path, gff_path):
    """Prepares, compresses, and indexes FASTA and GFF files specifically for VEP."""
    logging.info("Checking VEP reference index status...")
    
    vep_fna = f"{fna_path}.gz"
    vep_gff = f"{gff_path}.gz"

    # 1. Compress and Index FASTA for VEP
    if not os.path.exists(vep_fna) or not os.path.exists(f"{vep_fna}.fai"):
        run_cmd(f"bgzip -c {fna_path} > {vep_fna}", "Compressing FASTA for VEP")
        run_cmd(f"samtools faidx {vep_fna}", "Indexing compressed FASTA")

    # 2. Sort, Compress, and Index GFF3 for VEP
    if not os.path.exists(vep_gff) or not os.path.exists(f"{vep_gff}.tbi"):
        # VEP requires strict coordinate sorting and tabix indexing
        sort_cmd = (
            f"grep -v '#' {gff_path} | "
            f"sort -k1,1 -k4,4n -k5,5n -t$'\t' | "
            f"bgzip -c > {vep_gff}"
        )
        run_cmd(sort_cmd, "Coordinate-sorting and compressing GFF3 for VEP")
        run_cmd(f"tabix -p gff {vep_gff}", "Tabix indexing GFF3 file")

    return vep_fna, vep_gff

def parse_vep_vcf(vcf_path, csv_out_path):
    """Parses VEP's CSQ annotated VCF lines into a human-readable summary table."""
    if not os.path.exists(vcf_path):
        return

    with open(vcf_path, "r") as infile, open(csv_out_path, "w") as csvfile:
        csvfile.write("Position,Ref,Alt,Gene,Consequence,Impact\n")
        
        for line in infile:
            if line.startswith("#"):
                continue
            
            columns = line.strip().split("\t")
            pos = columns[1]
            ref = columns[3]
            alt = columns[4]
            info = columns[7]
            
            gene, consequence, impact = "Unknown", "Unknown", "Unknown"
            
            # Extract VEP's CSQ block
            if "CSQ=" in info:
                csq_segment = [x for x in info.split(";") if x.startswith("CSQ=")][0]
                first_transcript = csq_segment.replace("CSQ=", "").split(",")[0]
                csq_fields = first_transcript.split("|")
                
                # Standard VEP CSQ format: Allele|Consequence|IMPACT|SYMBOL|Gene...
                if len(csq_fields) >= 4:
                    consequence = csq_fields[1] if csq_fields[1] else "Unknown"
                    impact = csq_fields[2] if csq_fields[2] else "Unknown"
                    gene = csq_fields[3] if csq_fields[3] else "Intergenic"

            csvfile.write(f"{pos},{ref},{alt},{gene},{consequence},{impact}\n")

def process_illumina_pair(r1_file, r2_file, sample_name, pipeline_ref, vep_fna, vep_gff, threads, mem):
    """Handles quality trimming, alignment, variant calling, and VEP annotation."""
    logging.info(f"\n{'='*60}\n PROCESSING ILLUMINA SAMPLE: {sample_name}\n{'='*60}")

    # 1. Quality Trimming
    trim_r1 = f"alignment/{sample_name}_trimmed_R1.fq.gz"
    trim_r2 = f"alignment/{sample_name}_trimmed_R2.fq.gz"
    fastp_cmd = (
        f"fastp -i {r1_file} -I {r2_file} -o {trim_r1} -O {trim_r2} "
        f"--thread {threads} --html reports/{sample_name}_fastp.html "
        f"--detect_adapter_for_pe"
    )
    run_cmd(fastp_cmd, f"[{sample_name}] Fastp Quality Trimming")

    # 2. Alignment
    sorted_bam = f"alignment/{sample_name}_sorted.bam"
    align_cmd = (
        f"bwa mem -t {threads} {pipeline_ref} {trim_r1} {trim_r2} | "
        f"samtools view -@ {threads} -u - | "
        f"samtools sort -@ {threads} -m {mem} -o {sorted_bam} -"
    )
    run_cmd(align_cmd, f"[{sample_name}] BWA-MEM Alignment & Sorting")
    run_cmd(f"samtools index -@ {threads} {sorted_bam}", f"[{sample_name}] BAM Indexing")

    # 3. Variant Calling
    filtered_vcf = f"variants/{sample_name}_filtered.vcf"
    call_cmd = (
        f"bcftools mpileup --threads {threads} -Ou -f {pipeline_ref} {sorted_bam} | "
        f"bcftools call --threads {threads} -Ou -mv | "
        f"bcftools filter --threads {threads} -i 'QUAL>=20 && DP>=10' -o {filtered_vcf}"
    )
    run_cmd(call_cmd, f"[{sample_name}] Bcftools Variant Calling")

    # 4. Ensembl VEP Annotation
    final_vcf = f"results/batch_outputs/{sample_name}_vep_annotated.vcf"
    vep_cmd = (
        f"vep -i {filtered_vcf} -o {final_vcf} --vcf "
        f"--fasta {vep_fna} --gff {vep_gff} "
        f"--force_overwrite --everything --fork {threads}"
    )
    run_cmd(vep_cmd, f"[{sample_name}] VEP Functional Annotation")

    # 5. Parsing & Cleanup
    summary_csv = f"results/batch_outputs/{sample_name}_vep_summary.csv"
    parse_vep_vcf(final_vcf, summary_csv)

    logging.info(f"[{sample_name}] Cleaning up intermediate fastq files...")
    for f in [trim_r1, trim_r2]:
        if os.path.exists(f): os.remove(f)

# ==========================================
# EXECUTION LOGIC
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Universal Illumina VEP Batch Pipeline")
    parser.add_argument("--batch_dir", required=True, help="Directory containing paired fastq.gz files")
    parser.add_argument("--ref", required=True, help="Path to reference genome (.fna)")
    parser.add_argument("--gff", required=True, help="Path to matching reference annotation (.gff)")
    parser.add_argument("--threads", default="4", help="Number of CPU threads")
    parser.add_argument("--mem", default="1G", help="Memory allocation for samtools sort")
    
    args = parser.parse_args()

    for d in DIRS: 
        os.makedirs(d, exist_ok=True)
        
    setup_logging()
    logging.info("Starting automated Illumina paired-end VEP pipeline.")

    # --- SETUP PHASE ---
    pipeline_ref = f"reference/{os.path.basename(args.ref)}"
    if not os.path.exists(pipeline_ref): 
        shutil.copy(args.ref, pipeline_ref)

    if not os.path.exists(f"{pipeline_ref}.bwt"):
        run_cmd(f"bwa index {pipeline_ref}", "BWA Reference Indexing")

    # Prepare specific Bgzipped and Tabix-indexed files required by VEP
    vep_fna, vep_gff = prepare_vep_references(args.ref, args.gff)

    # --- PAIRED-END MATCHING PHASE ---
    r1_files = glob.glob(os.path.join(args.batch_dir, "*_1.fastq.gz")) + glob.glob(os.path.join(args.batch_dir, "*_R1.fastq.gz"))
    
    if not r1_files:
        logging.error(f"No R1 fastq files found in {args.batch_dir}. Ensure files end in _1.fastq.gz or _R1.fastq.gz")
        sys.exit(1)

    logging.info(f"Identified {len(r1_files)} paired samples for processing.")

    for r1 in r1_files:
        if "_1.fastq.gz" in r1:
            r2 = r1.replace("_1.fastq.gz", "_2.fastq.gz")
            sample_name = os.path.basename(r1).replace("_1.fastq.gz", "")
        else:
            r2 = r1.replace("_R1.fastq.gz", "_R2.fastq.gz")
            sample_name = os.path.basename(r1).replace("_R1.fastq.gz", "")

        if not os.path.exists(r2):
            logging.warning(f"Missing reverse pair (R2) for {r1}. Skipping {sample_name}.")
            continue
            
        process_illumina_pair(r1, r2, sample_name, pipeline_ref, vep_fna, vep_gff, args.threads, args.mem)
        
    logging.info("\n" + "="*60 + "\n ALL ILLUMINA SAMPLES PROCESSED WITH VEP \n" + "="*60)

if __name__ == "__main__":
    main()