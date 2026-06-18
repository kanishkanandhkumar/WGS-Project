#!/bin/bash

# Ensure a dedicated folder exists for our batch results
mkdir -p results/batch_outputs

echo "================================================================================"
echo " Starting Batch Viral Variant Pipeline Processing"
echo "================================================================================"

# Loop through every gzipped fastq file inside the data directory
for fastq_file in data/*.fastq.gz; do
    
    # Extract just the sample name (e.g., data/ERR16055133.fastq.gz -> ERR16055133)
    sample_name=$(basename "$fastq_file" .fastq.gz)
    
    echo ""
    echo "--------------------------------------------------------------------------------"
    echo " PROCESSING SAMPLE: $sample_name"
    echo "--------------------------------------------------------------------------------"
    
    # Run your established pipeline on the current sample
    python run_viral_amplicon_pipeline.py \
      --fastq "$fastq_file" \
      --ref reference/sars_cov_2.fna \
      --gff reference/sars_cov_2.gff
      
    # 1. Safely move and rename the CSV report
    if [ -f results/sars_cov_2_mutation_summary.csv ]; then
        mv results/sars_cov_2_mutation_summary.csv results/batch_outputs/${sample_name}_mutation_summary.csv
        echo "[SUCCESS] Saved CSV report: ${sample_name}_mutation_summary.csv"
    fi

    # 2. Safely move and rename ANY generated VCF files (including the SnpEff annotated ones)
    for vcf_file in results/*.vcf; do
        # Check if the file actually exists (prevents errors if no VCF is found)
        if [ -f "$vcf_file" ]; then
            vcf_basename=$(basename "$vcf_file")
            # Move and prefix the VCF with the sample name
            mv "$vcf_file" "results/batch_outputs/${sample_name}_${vcf_basename}"
            echo "[SUCCESS] Saved Annotated VCF: ${sample_name}_${vcf_basename}"
        fi
    done

done

echo ""
echo "================================================================================"
echo " All samples processed! Check results/batch_outputs/ for your VCF and CSV files."
echo "================================================================================"
