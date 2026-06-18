#!/usr/bin/env python3
import os
import sys
import subprocess
import logging
import time
import argparse

# ==========================================
# CONSTANTS & DIRECTORY CONFIGURATION
# ==========================================
DIRS = ["reference", "alignment", "variants", "results", "reports", "logs"]

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler("logs/nanopore_snpeff_pipeline.log"),
            logging.StreamHandler(sys.stdout)
        ]
    )

def run_cmd(cmd, step_name):
    start_time = time.time()
    logging.info(f"STARTING: {step_name}")
    try:
        subprocess.run(
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

def parse_snpeff_vcf(vcf_path, csv_out_path):
    """Parses SnpEff's ANN annotated VCF lines into a human-readable summary table."""
    if not os.path.exists(vcf_path):
        return

    with open(vcf_path, "r") as infile, open(csv_out_path, "w") as csvfile:
        csvfile.write("Position,Ref,Alt,Gene,Consequence,Impact,Amino_Acid_Change\n")
        
        for line in infile:
            if line.startswith("#"):
                continue
            
            columns = line.strip().split("\t")
            pos = columns[1]
            ref = columns[3]
            alt = columns[4]
            info = columns[7]
            
            gene, consequence, impact, aa_change = "Unknown", "Unknown", "Unknown", "Unknown"
            
            # Extract SnpEff's ANN block
            if "ANN=" in info:
                ann_segment = [x for x in info.split(";") if x.startswith("ANN=")][0]
                first_transcript = ann_segment.replace("ANN=", "").split(",")[0]
                ann_fields = first_transcript.split("|")
                
                # Standard SnpEff ANN format: Allele|Annotation|Annotation_Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank|HGVS.c|HGVS.p...
                if len(ann_fields) >= 11:
                    consequence = ann_fields[1] if ann_fields[1] else "Unknown"
                    impact = ann_fields[2] if ann_fields[2] else "Unknown"
                    gene = ann_fields[3] if ann_fields[3] else "Intergenic"
                    aa_change = ann_fields[10] if ann_fields[10] else "None"

            csvfile.write(f"{pos},{ref},{alt},{gene},{consequence},{impact},{aa_change}\n")

# ==========================================
# EXECUTION LOGIC
# ==========================================
def main():
    parser = argparse.ArgumentParser(description="Nanopore Long-Read Viral Pipeline (Minimap2 + SnpEff)")
    parser.add_argument("--fastq", required=True, help="Path to raw Nanopore fastq.gz file")
    parser.add_argument("--ref", required=True, help="Path to reference genome (.fna)")
    parser.add_argument("--gff", required=True, help="Path to reference annotation (.gff)")
    parser.add_argument("--threads", default="4", help="Number of CPU threads")
    
    args = parser.parse_args()

    for d in DIRS: 
        os.makedirs(d, exist_ok=True)
        
    setup_logging()
    logging.info("Starting automated viral Nanopore amplicon analysis pipeline.")
    
    sample_name = os.path.basename(args.fastq).replace(".fastq.gz", "").replace(".fq.gz", "")
    ref_basename = os.path.basename(args.ref)
    
    # 1. Reference Indexing
    if not os.path.exists(f"{args.ref}.fai"):
        run_cmd(f"samtools faidx {args.ref}", "Samtools FASTA indexing")
        
    mmi_index = f"reference/{ref_basename}.mmi"
    if not os.path.exists(mmi_index):
        run_cmd(f"minimap2 -d {mmi_index} {args.ref}", f"Generating Minimap2 index cache ({ref_basename}.mmi)")

    # 2. Nanopore Alignment (Minimap2)
    sorted_bam = f"alignment/{sample_name}_sorted.bam"
    align_cmd = (
        f"minimap2 -ax map-ont -t {args.threads} {mmi_index} {args.fastq} | "
        f"samtools view -@ {args.threads} -u - | "
        f"samtools sort -@ {args.threads} -o {sorted_bam} -"
    )
    run_cmd(align_cmd, "Minimap2 alignment and coordinate sorting")
    run_cmd(f"samtools index -@ {args.threads} {sorted_bam}", "BAM spatial file indexing")

    # 3. Variant Calling & Quality Filtering
    filtered_vcf = f"variants/{sample_name}_filtered.vcf"
    call_cmd = (
        f"bcftools mpileup --threads {args.threads} -Ou -f {args.ref} {sorted_bam} | "
        f"bcftools call --threads {args.threads} -Ou -mv | "
        f"bcftools filter --threads {args.threads} -i 'QUAL>=20 && DP>=10' -o {filtered_vcf}"
    )
    run_cmd(call_cmd, "Bcftools multi-threaded variant calling and quality filtering")

    # 4. Building Local SnpEff Database (Auto-builds if missing)
    snpeff_db_dir = "reference/snpeff_db"
    genome_version = ref_basename.replace('.fna', '')
    if not os.path.exists(f"{snpeff_db_dir}/{genome_version}/snpEffectPredictor.bin"):
        logging.info("Building local SnpEff database for the custom reference...")
        os.makedirs(f"{snpeff_db_dir}/{genome_version}", exist_ok=True)
        # Copy references into the structure SnpEff expects
        subprocess.run(f"cp {args.ref} {snpeff_db_dir}/{genome_version}/sequences.fa", shell=True)
        subprocess.run(f"cp {args.gff} {snpeff_db_dir}/{genome_version}/genes.gff", shell=True)
        # Create a temporary config
        with open("snpeff.config", "w") as f:
            f.write(f"data.dir = {os.path.abspath(snpeff_db_dir)}\n")
            f.write(f"{genome_version}.genome : Custom Viral Genome\n")
        # Build the database
        run_cmd(f"snpEff build -c snpeff.config -gff3 -v -noCheckCds -noCheckProtein {genome_version}", "Compiling custom SnpEff genomic database")

    # 5. SnpEff Annotation
    annotated_vcf = f"results/{sample_name}_snpeff_annotated.vcf"
    snpeff_cmd = (
        f"snpEff -c snpeff.config {genome_version} {filtered_vcf} > {annotated_vcf}"
    )
    run_cmd(snpeff_cmd, "SnpEff Functional Annotation")

    # 6. Parse Output to CSV
    summary_csv = f"results/{sample_name}_mutation_summary.csv"
    parse_snpeff_vcf(annotated_vcf, summary_csv)
    
    logging.info(f"Pipeline complete for {sample_name}.")

if __name__ == "__main__":
    main()