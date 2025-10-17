#!/usr/bin/env python3
"""
Main script for analyzing chromatin loop lengths from BEDPE files.

This script reads a configuration YAML file and performs loop length analysis,
motif scanning, and visualization of the results.
"""

import os
import sys
import pandas as pd
from src.config import load_config, get_motif_config, get_apa_config
from src.utils import (
    read_bedpe,
    filter_intrachromosomal,
    filter_interchromosomal,
    filter_by_loop_length,
    setup_output_directories
)
from src.visualizations import plot_loop_length_violin_by_category, plot_apa_comparison_heatmaps
from src.motif_analysis import (
    extract_anchors_to_bed,
    extract_anchor_fasta,
    run_fimo,
    parse_fimo_results,
    determine_anchor_states,
    classify_loops,
    generate_statistics,
    output_separated_bedpe
)
from src.apa_analysis import run_apa_by_category


def main():
    """
    Main function to orchestrate the loop length analysis workflow.
    """
    # Check for config file argument
    if len(sys.argv) != 2:
        print("[main][error] Usage: python main.py <config.yaml>")
        print("[main][error] See configs/example_config.yaml for configuration template")
        sys.exit(1)

    config_path = sys.argv[1]

    # Load configuration
    print("[main][info] Loading configuration...")
    config = load_config(config_path)

    # Extract config values
    bedpe_path = config['input']['bedpe']
    has_header = config['input'].get('has_header', True)
    output_dir = config['output']['directory']
    include_interchromosomal = config['loops'].get('include_interchromosomal', False)

    print("[main][info] Starting loop length analysis...")
    print(f"[main][info] Input BEDPE file: {bedpe_path}")
    print(f"[main][info] Output directory: {output_dir}")

    # Setup output directory structure
    out_dirs = setup_output_directories(output_dir)

    # Read BEDPE file
    print("[main][info] Reading BEDPE file...")
    df = read_bedpe(bedpe_path, has_header=has_header)

    # Filter loops based on chromosome type
    if include_interchromosomal:
        print("[main][info] Using all loops (intra- and inter-chromosomal)")
        loops_df = df.copy()
    else:
        print("[main][info] Filtering to intra-chromosomal loops only")
        loops_df = filter_intrachromosomal(df)

    # Motif analysis: extract anchors and FASTA if requested
    motif_config = get_motif_config(config)
    if motif_config:
        print("[main][info] Starting motif analysis workflow...")

        # Extract anchors to BED file
        print("[main][info] Extracting anchors to BED format...")
        bed_file = os.path.join(out_dirs['temp'], 'anchors.bed')
        extract_anchors_to_bed(loops_df, bed_file)

        # Extract anchor sequences to FASTA
        print("[main][info] Extracting anchor sequences to FASTA...")
        fasta_file = os.path.join(out_dirs['temp'], 'anchors.fa')
        success = extract_anchor_fasta(bed_file, motif_config['fasta'], fasta_file)

        if not success:
            print("[main][error] Failed to extract anchor FASTA sequences")
            sys.exit(1)

        # Run FIMO motif scanning
        print("[main][info] Running FIMO motif scanning...")
        fimo_output_dir = os.path.join(out_dirs['temp'], 'fimo')
        os.makedirs(fimo_output_dir, exist_ok=True)

        fimo_results = run_fimo(
            fasta_file,
            motif_config['motif'],
            fimo_output_dir,
            fimo_threshold=motif_config['fimo_threshold']
        )

        if not fimo_results:
            print("[main][error] FIMO motif scanning failed or found no hits")
            sys.exit(1)

        # Parse FIMO results
        print("[main][info] Parsing FIMO results...")
        anchor_hits = parse_fimo_results(fimo_results)

        # Determine anchor states
        print("[main][info] Determining anchor states...")
        anchor_states = determine_anchor_states(
            anchor_hits,
            loops_df,
            fdr_threshold=motif_config['fdr_threshold'],
            optimize_convergent=motif_config['optimize_convergent']
        )

        # Classify loops
        print("[main][info] Classifying loops...")
        classification_result = classify_loops(loops_df, anchor_states)
        classifications = classification_result['classifications']
        category_counts = classification_result['category_counts']

        # Generate statistics
        print("[main][info] Generating classification statistics...")
        stats = generate_statistics(loops_df, classifications, category_counts, out_dirs['generated_files'])

        # Output separated BEDPE files
        print("[main][info] Writing separated BEDPE files...")
        bedpe_files = output_separated_bedpe(loops_df, classifications, anchor_states, out_dirs['generated_files'])

        # Generate violin plot
        print("[main][info] Generating violin plot by category...")

        # Load all category BEDPE files
        category_files = {
            'all_loops': loops_df,
            'zero_anchors': os.path.join(out_dirs['generated_files'], 'zero_anchors.bedpe'),
            'one_anchor': os.path.join(out_dirs['generated_files'], 'one_anchor.bedpe'),
            'all_two_anchors': os.path.join(out_dirs['generated_files'], 'all_two_anchors.bedpe'),
            'convergent': os.path.join(out_dirs['generated_files'], 'convergent.bedpe'),
            'tandem_forward': os.path.join(out_dirs['generated_files'], 'tandem_forward.bedpe'),
            'tandem_reverse': os.path.join(out_dirs['generated_files'], 'tandem_reverse.bedpe'),
            'divergent': os.path.join(out_dirs['generated_files'], 'divergent.bedpe')
        }

        # Read all category dataframes
        bedpe_dfs = {}
        for category, path_or_df in category_files.items():
            if isinstance(path_or_df, str):
                # It's a file path
                if os.path.exists(path_or_df):
                    bedpe_dfs[category] = pd.read_csv(path_or_df, sep='\t')
                else:
                    bedpe_dfs[category] = pd.DataFrame()  # Empty dataframe if file doesn't exist
            else:
                # It's already a dataframe
                bedpe_dfs[category] = path_or_df

        # Generate violin plot
        plot_loop_length_violin_by_category(
            bedpe_dfs,
            out_dirs['figures'],
            generated_files_dir=out_dirs['generated_files'],
            filename='loop_length_violin_by_category.png'
        )
    
    # APA analysis if configured
    apa_config = get_apa_config(config)
    if apa_config:
        print("[main][info] Starting APA (Aggregate Peak Analysis) workflow...")

        # Filter loops by minimum length for APA if specified
        apa_bedpe_dfs = bedpe_dfs.copy()
        min_loop_length = apa_config.get('min_loop_length')
        if min_loop_length is not None:
            print(f"[main][info] Filtering loops for APA (min length: {min_loop_length/1000:.0f}kb)")
            apa_bedpe_dfs = {}
            for category, loops_df in bedpe_dfs.items():
                if len(loops_df) == 0:
                    apa_bedpe_dfs[category] = loops_df
                    remaining = 0
                else:
                    filtered_df = filter_by_loop_length(loops_df, min_loop_length)
                    apa_bedpe_dfs[category] = filtered_df
                    remaining = len(filtered_df)

                original = len(loops_df)
                removed = original - remaining
                print(f"[main][info]   {category:20s} {original:5d} → {remaining:5d} loops ({removed} removed)")

        # Copy Hi-C file to temp before processing (do this once, not per category)
        import shutil
        hic_file = apa_config['hic_file']
        hic_temp_file = os.path.join(out_dirs['temp'], 'hic_apa.cool')
        if not os.path.exists(hic_temp_file):
            print(f"[main][info] Copying Hi-C file to temp: {hic_temp_file}")
            shutil.copy2(hic_file, hic_temp_file)
            # Update config to use temp copy
            apa_config['hic_file'] = hic_temp_file

        # Run APA for each loop category
        apa_results = run_apa_by_category(
            apa_bedpe_dfs,
            apa_config,
            out_dirs['generated_files']
        )

        # Generate comparison visualization
        print("[main][info] Generating APA comparison heatmaps...")
        plot_apa_comparison_heatmaps(
            apa_results,
            out_dirs['figures'],
            resolution=apa_config['resolution'],
            window_size=apa_config['window_size']
        )

        print("[main][info] APA analysis complete")

    # Print summary
    print("\n[main][info] Analysis Complete!")
    print(f"[main][info] Output files saved to: {output_dir}")
    print(f"[main][info]   Figures: {out_dirs['figures']}")
    print(f"[main][info]   Generated files: {out_dirs['generated_files']}")
    print(f"[main][info]   Temporary files: {out_dirs['temp']}")


if __name__ == '__main__':
    main()
