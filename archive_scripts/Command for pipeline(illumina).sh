# Auto-detect your reference genome files (.fna and .gff)
REF_FNA=$(ls reference/*.fna | head -n 1)
REF_GFF=$(ls reference/*.gff | head -n 1)

# Run the Illumina pipeline
python run_illumina_batch_pipeline.py \
    --batch_dir data/ \
    --ref "$REF_FNA" \
    --gff "$REF_GFF" \
    --threads 4