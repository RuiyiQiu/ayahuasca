import sys
import os
import argparse
import numpy as np
import pandas as pd

def find_precise_boundaries_with_prior(chrom_df, chrom, prior_start, prior_end, 
                                       window_size=50000, step_size=5000, threshold=0.5, consecutive_windows=1):
    """
    Advanced sliding window algorithm that is resilient to internal low-signal valleys/gaps.
    Now supports fully dynamic windows, steps, and consecutive requirements via CLI.
    """
    # Safety buffer zone extension (1Mb on both sides)
    search_start = max(0, prior_start - 1000000)
    search_end = prior_end + 1000000
    
    region_df = chrom_df[(chrom_df['end'] >= search_start) & (chrom_df['start'] <= search_end)].copy()
    
    if region_df.empty:
        print(f"  [Warning] No signal data found for {chrom}. Skipping micro-trimming.")
        return prior_start, prior_end

    # Build high-resolution sliding grid positions
    grid_positions = np.arange(search_start, search_end, step_size)
    grid_values = []
    
    # Calculate weighted mean signal for each window
    for pos in grid_positions:
        w_start = pos
        w_end = pos + window_size
        overlap = region_df[(region_df['end'] > w_start) & (region_df['start'] < w_end)]
        if not overlap.empty:
            overlap_start = np.maximum(overlap['start'], w_start)
            overlap_end = np.minimum(overlap['end'], w_end)
            lengths = overlap_end - overlap_start
            weighted_mean = np.sum(overlap['value'] * lengths) / np.sum(lengths)
            grid_values.append(weighted_mean)
        else:
            grid_values.append(0.0)
            
    grid_values = np.array(grid_values)
    
    # 1. Locate precise Start Boundary (First rise above threshold)
    is_above_threshold = grid_values >= threshold
    kernel = np.ones(consecutive_windows)
    consecutive_counts = np.convolve(is_above_threshold.astype(int), kernel, mode='valid')
    
    valid_starts = np.where(consecutive_counts == consecutive_windows)[0]
    if len(valid_starts) == 0:
        print(f"  [Warning] Cutoff {threshold} never met consecutively ({consecutive_windows} windows) on {chrom}. Reverting to MACS3.")
        return prior_start, prior_end
        
    start_idx = valid_starts[0]
    precise_start = int(grid_positions[start_idx])
    
    # 2. NEW FIX LOGIC FOR END BOUNDARY:
    valid_ends = np.where(consecutive_counts == consecutive_windows)[0]
    
    # The last index where N consecutive windows were ABOVE threshold
    last_high_idx = valid_ends[-1] if len(valid_ends) > 0 else start_idx
    
    # Add the consecutive window offset to get the end position of the last active window
    precise_end = int(grid_positions[last_high_idx] + window_size)
    
    # Safety boundary capping based on prior range
    if precise_end > search_end:
        precise_end = prior_end
        
    print(f"  -> Refined {chrom}: {prior_start:,}-{prior_end:,} ===> {precise_start:,}-{precise_end:,} bp")
    return precise_start, precise_end


if __name__ == "__main__":
    # Setup argument parser to match standard bash execution
    parser = argparse.ArgumentParser(description="Micro-trim centromere boundaries using sliding window method.")
    
    # Position arguments (Required)
    parser.add_argument("-i", "--input_bed", help="Path to the filtered MACS3 input file (MACS3.flt.out) in BED format.")
    parser.add_argument("-b", "--bedgraph", 
                        help="Path to the original deeptools bamCompare signal file. (Default: compare.bedGraph)")
    parser.add_argument("-o", "--output_file", help="Path to save the final refined output boundaries (cen_precise_locate.out).")
    
    # Optional arguments (With advanced dynamic adjustments)
    parser.add_argument("-t", "--threshold", type=float, default=0.5, 
                        help="Log2 ChIP/Input threshold value for sliding window. (Default: 0.5)")
    parser.add_argument("-w", "--window_size", type=int, default=10000,
                        help="Sliding window size in bp. (Default: 10000)")
    parser.add_argument("-s", "--step_size", type=int, default=2000,
                        help="Sliding step size in bp. (Default: 2000)")
    parser.add_argument("-c", "--consecutive", type=int, default=5,
                        help="Number of consecutive windows required to meet threshold. (Default: 5)")
    
    args = parser.parse_args()
    
    # Validate the existence of required input genomic files
    if not os.path.exists(args.input_bed):
        print(f"Error: The input BED file '{args.input_bed}' does not exist.")
        sys.exit(1)
    if not os.path.exists(args.bedgraph):
        print(f"Error: The signal BedGraph file '{args.bedgraph}' does not exist. Please specify it using '-b <file>'.")
        sys.exit(1)
        
    print(f"Step 1: Loading continuous signal data from {args.bedgraph} ...")
    signal_df = pd.read_csv(args.bedgraph, sep='\t', header=None, names=['chrom', 'start', 'end', 'value'])
    
    print(f"Step 2: Loading macro prior boundaries from {args.input_bed} ...")
    prior_df = pd.read_csv(args.input_bed, sep='\t', header=None, names=['chrom', 'start', 'end'])
    
    final_refined_rows = []
    
    print(f"\nStep 3: Sliding-window micro-trimming (Window: {args.window_size}bp, Step: {args.step_size}bp, Consecutive: {args.consecutive}, Threshold: {args.threshold})...")
    
    for _, row in prior_df.iterrows():
        ch = row['chrom']
        p_start = int(row['start'])
        p_end = int(row['end'])
        
        chrom_signal = signal_df[signal_df['chrom'] == ch]
        
        # Execute the optimization locator with dynamic parameters from user CLI inputs
        ref_start, ref_end = find_precise_boundaries_with_prior(
            chrom_df=chrom_signal,
            chrom=ch,
            prior_start=p_start,
            prior_end=p_end,
            window_size=args.window_size,
            step_size=args.step_size,
            threshold=args.threshold,
            consecutive_windows=args.consecutive
        )
        
        final_refined_rows.append([ch, ref_start, ref_end])
        
    # Write the high-precision results back into the customized output file path
    refined_df = pd.DataFrame(final_refined_rows)
    refined_df.to_csv(args.output_file, sep='\t', header=False, index=False)
    
    print(f"\n[Success] High-precision centromere localization completed!")
    print(f"Final coordinates are successfully saved to: '{args.output_file}'")