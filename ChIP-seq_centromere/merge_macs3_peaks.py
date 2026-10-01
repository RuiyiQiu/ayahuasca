import sys
import os
import argparse
import pandas as pd

def merge_dense_peaks_all_chroms(macs_path, max_gap=100000):
    
    # Read MACS3 output file (compatible with both narrowPeak and broadPeak)
    # Only extract the first 3 columns: chrom, start, end
    df = pd.read_csv(macs_path, sep='\t', header=None, comment='#')
    df = df[[0, 1, 2]]
    df.columns = ['chrom', 'start', 'end']
    
    # Get a sorted list of all unique chromosomes present in the file
    unique_chroms = sorted(df['chrom'].unique())
    print(f"Found {len(unique_chroms)} unique chromosomes/scaffolds to process.\n")
    
    # Initialize a list to store final centromere boundaries for each chromosome
    final_results = []
    
    # Loop through every single chromosome automatically
    for chrom in unique_chroms:
        # Filter and sort rows belonging to the current chromosome
        chrom_df = df[df['chrom'] == chrom].sort_values(by='start').copy()
        
        merged_peaks = []
        current_start = None
        current_end = None
        
        # Core agglomerative clustering logic based on genomic distance
        for _, row in chrom_df.iterrows():
            if current_start is None:
                current_start = row['start']
                current_end = row['end']
            else:
                # If the gap between current peak start and previous peak end <= max_gap, merge them
                if row['start'] - current_end <= max_gap:
                    current_end = max(current_end, row['end'])
                else:
                    # Save the previous completed cluster region
                    merged_peaks.append((current_start, current_end))
                    current_start = row['start']
                    current_end = row['end']
                    
        # Append the final remaining cluster for this chromosome
        if current_start is not None:
            merged_peaks.append((current_start, current_end))
        
        # Find the longest merged region (the primary centromere candidate cluster)
        if merged_peaks:
            merged_df = pd.DataFrame(merged_peaks, columns=['start', 'end'])
            merged_df['length'] = merged_df['end'] - merged_df['start']
            
            # Sort by length descending and pick the top 1 largest block
            cen_cluster = merged_df.sort_values(by='length', ascending=False).iloc[0]
            
            final_results.append({
                'chrom': chrom,
                'macs_start': int(cen_cluster['start']),
                'macs_end': int(cen_cluster['end']),
                'total_span_kb': abs(int(cen_cluster['end']) - int(cen_cluster['start'])) / 1000.0
            })
            
    # Convert all gathered results into a structured pandas DataFrame
    results_df = pd.DataFrame(final_results)
    return results_df

if __name__ == "__main__":
    # Setup command-line argument parser for high flexibility
    parser = argparse.ArgumentParser(description="Batch merge MACS3 dense peaks and define centromere boundaries.")
    
    # Position argument: Input file path
    parser.add_argument("input_file", help="Path to the MACS3 broadPeak or narrowPeak file.")
    
    # Optional argument: Output BED file path (Defaults to 'macs3_merged_centromeres.bed')
    parser.add_argument("-o", "--output", default="macs3_merged_centromeres.bed", 
                        help="Path to the output BED file. (Default: macs3_merged_centromeres.bed)")
    
    # Optional argument: Maximum allowed gap to merge two peaks (Defaults to 100000 bp)
    parser.add_argument("-g", "--max_gap", type=int, default=100000, 
                        help="Maximum genomic distance (bp) allowed between peaks to merge them. (Default: 100000)")
    
    # Parse all arguments from user input
    args = parser.parse_args()
    
    # Validate the existence of input file
    if not os.path.exists(args.input_file):
        print(f"Error: The input file '{args.input_file}' does not exist.")
        sys.exit(1)
        
    # Run the batch merging pipeline with specified max gap
    summary_df = merge_dense_peaks_all_chroms(args.input_file, max_gap=args.max_gap)
    
    # Save the results into the custom defined output BED path
    summary_df[['chrom', 'macs_start', 'macs_end']].to_csv(args.output, sep='\t', header=False, index=False)
    print(f"\n[Success] Batch results successfully saved to custom BED file: '{args.output}'")
