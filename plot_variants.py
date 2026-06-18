import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import glob

def main():
    parser = argparse.ArgumentParser(description="Plot variant metrics from WGS batch outputs.")
    parser.add_argument("--dir", required=True, help="Directory containing mutation summary CSVs")
    args = parser.parse_args()

    print(f"[*] Scanning target directory: {args.dir}")
    csv_files = glob.glob(os.path.join(args.dir, "*_mutation_summary.csv"))
    
    if not csv_files:
        print(f"[!] Error: No mutation summary CSV files found in directory '{args.dir}'")
        return

    for csv_file in csv_files:
        print(f"[*] Loading dataset: {csv_file}")
        try:
            df = pd.read_csv(csv_file)
            if df.empty:
                print(f"[!] Warning: {csv_file} is empty. Skipping.")
                continue
            
            print(f"[*] Detected columns: {list(df.columns)}")
            
            plt.figure(figsize=(10, 6))
            
            # Dynamically look for common variant annotation columns
            target_col = None
            potential_cols = ['Impact', 'Effect', 'Type', 'Mutation_Type', 'Functional_Class', 'Chromosome', 'Chr']
            for col in potential_cols:
                if col in df.columns:
                    target_col = col
                    break
            
            if not target_col:
                text_cols = [c for c in df.columns if df[c].dtype == 'object']
                if text_cols:
                    target_col = text_cols[0]
            
            if target_col:
                print(f"[*] Generating metric distribution chart for column: '{target_col}'")
                counts = df[target_col].value_counts()
                counts.plot(kind='bar', color='#2b5c8f', edgecolor='black', alpha=0.85)
                plt.title(f"Variant Distribution Analysis ({target_col})\nDataset: {os.path.basename(csv_file)}", fontsize=12, fontweight='bold')
                plt.xlabel(target_col, fontsize=10)
                plt.ylabel("Variant Count", fontsize=10)
                plt.grid(axis='y', linestyle='--', alpha=0.5)
                plt.xticks(rotation=45, ha='right')
            else:
                print("[*] No categorical data found. Plotting raw metric matrix rows.")
                df.plot(kind='line', marker='o')
                plt.title(f"Genomic Variant Matrix Trends\nDataset: {os.path.basename(csv_file)}")
            
            plt.tight_layout()
            output_plot = os.path.join(args.dir, "batch_variant_summary_plot.png")
            plt.savefig(output_plot, dpi=300)
            plt.close()
            print(f"[+] Success! Visualization saved to: {output_plot}")
            
        except Exception as e:
            print(f"[!] Error processing data: {e}")

if __name__ == "__main__":
    main()
