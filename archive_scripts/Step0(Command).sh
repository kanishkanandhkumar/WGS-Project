# 1. Navigate to the project folder and activate the environment
conda activate wgs_pipeline
cd ~/Desktop/wgs_project

# 2. Safely wipe old workspace tracks (clears alignments, variants, reports, and results)
rm -rf alignment/* variants/* reports/* results/*
mkdir -p results/batch_outputs