# Auto-detect your reference genome files (.fna and .gff)
REF_FNA=$(ls reference/*.fna | head -n 1)
REF_GFF=$(ls reference/*.gff | head -n 1)

# Iterate through every fastq.gz file in the data folder
for fastq_file in data/*.fastq.gz; do
    
    # Isolate sample prefix name
    sample_name=$(basename "$fastq_file" .fastq.gz)
    
    echo "============================================================"
    echo " PROCESSING NANOPORE SAMPLE: $sample_name"
    echo "============================================================"
    
    # Execute the Nanopore pipeline
    python run_viral_amplicon_pipeline.py \
      --fastq "$fastq_file" \
      --ref "$REF_FNA" \
      --gff "$REF_GFF" \
      --threads 4
      
    # Move and isolate the mutation CSV report
    if [ -f results/*_mutation_summary.csv ]; then
        mv results/*_mutation_summary.csv results/batch_outputs/${sample_name}_mutation_summary.csv
        echo "[SAVED] Spreadsheet summary for $sample_name"
    fi

    # Move and isolate the final SnpEff annotated VCF file
    for vcf_file in results/*.vcf; do
        if [ -f "$vcf_file" ]; then
            vcf_basename=$(basename "$vcf_file")
            mv "$vcf_file" "results/batch_outputs/${sample_name}_${vcf_basename}"
            echo "[SAVED] Annotated VCF for $sample_name"
        fi
    done

done