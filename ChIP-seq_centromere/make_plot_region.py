import sys
import os
import argparse
import pandas as pd

def generate_plotting_regions(input_bed, output_bed):
    """
    Reads refined centromere coordinates, calculates length-dependent flanking regions 
    (1/4 of centromere length on each side), and outputs an optimized BED file for pyGenomeTracks.
    """
    print(f"Reading precise centromere coordinates from: {input_bed}")
    
    # Load the high-precision boundaries generated from the second script
    df = pd.read_csv(input_bed, sep='\t', header=None, names=['chrom', 'start', 'end'])
    
    plotting_rows = []
    
    # Process each chromosome dynamically
    for _, row in df.iterrows():
        chrom = row['chrom']
        cen_start = int(row['start'])
        cen_end = int(row['end'])
        
        # Calculate the exact core length of this specific centromere
        cen_length = cen_end - cen_start
        
        # Calculate flanking buffer size (exactly 1/4 of the core length)
        # Using integer division to keep coordinate values clean
        flank_size = cen_length // 4
        
        # Derive new absolute coordinates for full-range plotting visualization
        plot_start = max(0, cen_start - flank_size)
        plot_end = cen_end + flank_size
        
        print(f" -> {chrom}: Core Length = {cen_length/1000:,.1f} kb | Flank Size = {flank_size/1000:,.1f} kb")
        print(f"    Plotting Window: {chrom}:{plot_start:,}-{plot_end:,} bp")
        
        plotting_rows.append([chrom, plot_start, plot_end])
        
    # Convert results into a standard BED4 format DataFrame
    plot_df = pd.DataFrame(plotting_rows)
    plot_df.to_csv(output_bed, sep='\t', header=False, index=False)
    
    print(f"\n[Success] PyGenomeTracks visualization regions saved to: '{output_bed}'")


if __name__ == "__main__":
    # Setup command-line user interface
    parser = argparse.ArgumentParser(description="Create custom flank-extended BED files for high-impact pyGenomeTracks plots.")
    
    # Position arguments for standard pipeline running
    parser.add_argument("input_file", help="Path to the precise centromere BED file (e.g., cen_precise_locate.out).")
    parser.add_argument("output_file", help="Path to save the generated pyGenomeTracks region BED file.")
    
    args = parser.parse_args()
    
    # File execution safety check
    if not os.path.exists(args.input_file):
        print(f"Error: Input file '{args.input_file}' not found.")
        sys.exit(1)
        
    # Execute the processing function
    generate_plotting_regions(args.input_file, args.output_file)