"""
Motif analysis functions for CTCF loop analysis.

Handles anchor extraction, FASTA generation, FIMO motif scanning, and result parsing.
"""

import pandas as pd
import subprocess
import os


def extract_anchors_to_bed(bedpe_df, output_bed):
    """
    Extract loop anchors from BEDPE DataFrame and write to BED format.

    Each loop has two anchors (anchor1 and anchor2), which are written as separate BED records.
    Normalizes chromosome names by adding 'chr' prefix if missing.

    Parameters:
    -----------
    bedpe_df : pd.DataFrame
        DataFrame with columns: chr1, start1, end1, chr2, start2, end2
    output_bed : str
        Path to output BED file

    Returns:
    --------
    None
        Writes BED file to output_bed path
    """
    bed_records = []

    for loop_id, (idx, row) in enumerate(bedpe_df.iterrows(), start=1):
        # Use loop_id sequentially starting from 1, not the DataFrame index

        # Normalize chromosome names (add 'chr' prefix if missing)
        chr1 = str(row['chr1'])
        chr2 = str(row['chr2'])

        if not chr1.startswith('chr'):
            chr1 = f'chr{chr1}'
        if not chr2.startswith('chr'):
            chr2 = f'chr{chr2}'

        # Anchor 1
        anchor1_record = {
            'chrom': chr1,
            'start': int(row['start1']),
            'end': int(row['end1']),
            'name': f'L{loop_id}A1'
        }
        bed_records.append(anchor1_record)

        # Anchor 2
        anchor2_record = {
            'chrom': chr2,
            'start': int(row['start2']),
            'end': int(row['end2']),
            'name': f'L{loop_id}A2'
        }
        bed_records.append(anchor2_record)

    # Create BED DataFrame and write
    bed_df = pd.DataFrame(bed_records)
    bed_df.to_csv(output_bed, sep='\t', header=False, index=False)

    num_anchors = len(bed_df)
    print(f"[motif_analysis][info] Extracted {num_anchors} anchors to BED file: {output_bed}")


def extract_anchor_fasta(bed_file, fasta_ref, output_fasta):
    """
    Extract FASTA sequences for anchor regions using bedtools getfasta.

    Parameters:
    -----------
    bed_file : str
        Path to BED file with anchor coordinates
    fasta_ref : str
        Path to reference genome FASTA file
    output_fasta : str
        Path to output FASTA file

    Returns:
    --------
    bool
        True if successful, False if bedtools command failed
    """
    cmd = ['bedtools', 'getfasta', '-fi', fasta_ref, '-bed', bed_file, '-fo', output_fasta, '-nameOnly']

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"[motif_analysis][info] Extracted anchor sequences to FASTA: {output_fasta}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"[motif_analysis][error] bedtools getfasta failed: {e.stderr}")
        return False
    except FileNotFoundError:
        print(f"[motif_analysis][error] bedtools not found in PATH")
        return False


def run_fimo(fasta_file, motif_file, output_dir, fimo_threshold='1e-4'):
    """
    Run FIMO (Find Individual Motif Occurrences) for motif scanning.

    Parameters:
    -----------
    fasta_file : str
        Path to FASTA file with sequences to scan
    motif_file : str
        Path to motif file in MEME format
    output_dir : str
        Output directory for FIMO results
    fimo_threshold : str, default='1e-4'
        FIMO p-value threshold (e.g., '1e-4', '0.05')

    Returns:
    --------
    str or None
        Path to fimo.tsv results file if successful, None if failed
    """
    cmd = ['fimo', '--oc', output_dir, '--thresh', fimo_threshold, motif_file, fasta_file]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        fimo_results = os.path.join(output_dir, 'fimo.tsv')

        # Check if results file exists and has content
        if os.path.exists(fimo_results):
            with open(fimo_results, 'r') as f:
                lines = f.readlines()
            num_hits = len(lines) - 1  # Subtract header line

            if num_hits > 0:
                print(f"[motif_analysis][info] FIMO completed. Found {num_hits} motif hits")
                print(f"[motif_analysis][info] Results saved to: {fimo_results}")
                return fimo_results
            else:
                print(f"[motif_analysis][warning] FIMO completed but found no motif hits")
                return None
        else:
            print(f"[motif_analysis][error] FIMO output file not found: {fimo_results}")
            return None

    except subprocess.CalledProcessError as e:
        print(f"[motif_analysis][error] FIMO failed: {e.stderr}")
        return None
    except FileNotFoundError:
        print(f"[motif_analysis][error] fimo not found in PATH")
        return None


def parse_fimo_results_all_hits(fimo_tsv_path, fdr_threshold=0.05):
    """
    Parse FIMO results and extract ALL hits for each anchor.

    Parameters:
    -----------
    fimo_tsv_path : str
        Path to fimo.tsv results file
    fdr_threshold : float, default=0.05
        P-value threshold for filtering hits

    Returns:
    --------
    dict
        Dictionary mapping anchor IDs to list of all hits:
        {
            'L1A1': [
                {'start': int, 'end': int, 'pvalue': float, 'strand': '+', 'qvalue': float},
                {'start': int, 'end': int, 'pvalue': float, 'strand': '-', 'qvalue': float},
                ...
            ],
            ...
        }
    """
    anchor_all_hits = {}

    if not os.path.exists(fimo_tsv_path):
        print(f"[motif_analysis][error] FIMO results file not found: {fimo_tsv_path}")
        return anchor_all_hits

    # Read FIMO results with proper header handling, skip comment lines
    try:
        fimo_df = pd.read_csv(fimo_tsv_path, sep='\t', comment='#')
    except Exception as e:
        print(f"[motif_analysis][error] Failed to read FIMO TSV file: {e}")
        return {}

    # Verify we have the required columns
    required_cols = ['sequence_name', 'start', 'stop', 'strand', 'p-value']
    missing_cols = [col for col in required_cols if col not in fimo_df.columns]
    if missing_cols:
        print(f"[motif_analysis][error] Missing columns in FIMO output: {missing_cols}")
        print(f"[motif_analysis][error] Available columns: {list(fimo_df.columns)}")
        return {}

    # Check for NaN values in required columns BEFORE processing
    nan_mask = fimo_df[required_cols].isnull().any(axis=1)
    if nan_mask.any():
        num_nan_rows = nan_mask.sum()
        print(f"[motif_analysis][warning] Found {num_nan_rows} rows with NaN values in required columns")
        print(f"[motif_analysis][warning] Removing these rows before parsing")
        fimo_df = fimo_df.dropna(subset=required_cols)

    if len(fimo_df) == 0:
        print(f"[motif_analysis][error] No valid FIMO results after filtering out NaN rows")
        return {}

    # Parse each FIMO hit and store ALL hits (filtered by threshold)
    try:
        for _, row in fimo_df.iterrows():
            anchor_id = row['sequence_name']
            motif_start = int(row['start'])
            motif_end = int(row['stop'])
            strand = row['strand']
            pvalue = float(row['p-value'])
            qvalue = float(row['q-value']) if 'q-value' in row else pvalue

            # Filter by threshold
            if pvalue > fdr_threshold:
                continue

            # Initialize anchor entry if not exists
            if anchor_id not in anchor_all_hits:
                anchor_all_hits[anchor_id] = []

            # Store this hit
            anchor_all_hits[anchor_id].append({
                'start': motif_start,
                'end': motif_end,
                'pvalue': pvalue,
                'qvalue': qvalue,
                'strand': strand
            })

    except (ValueError, TypeError) as e:
        print(f"[motif_analysis][error] Failed to parse FIMO row values: {e}")
        print(f"[motif_analysis][error] Check that start, stop, and p-value columns contain numeric values")
        return {}

    num_anchors_with_hits = len(anchor_all_hits)
    total_hits = sum(len(hits) for hits in anchor_all_hits.values())
    print(f"[motif_analysis][info] Parsed ALL FIMO hits: {num_anchors_with_hits} anchors with hits, {total_hits} total hits")

    return anchor_all_hits


def parse_fimo_results(fimo_tsv_path):
    """
    Parse FIMO results and extract best forward and reverse hits for each anchor.

    For each anchor, retains at most TWO hits:
    - Best forward strand hit (lowest p-value on + strand)
    - Best reverse strand hit (lowest p-value on - strand)

    Parameters:
    -----------
    fimo_tsv_path : str
        Path to fimo.tsv results file

    Returns:
    --------
    dict
        Dictionary mapping anchor IDs to their best hits:
        {
            'L1A1': {
                'forward': {'start': int, 'end': int, 'pvalue': float, 'strand': '+'},
                'reverse': {'start': int, 'end': int, 'pvalue': float, 'strand': '-'}
            },
            ...
        }
        Where forward/reverse can be None if no hits found for that strand.
    """
    anchor_hits = {}

    if not os.path.exists(fimo_tsv_path):
        print(f"[motif_analysis][error] FIMO results file not found: {fimo_tsv_path}")
        return anchor_hits

    # Read FIMO results with proper header handling, skip comment lines
    try:
        fimo_df = pd.read_csv(fimo_tsv_path, sep='\t', comment='#')
    except Exception as e:
        print(f"[motif_analysis][error] Failed to read FIMO TSV file: {e}")
        return {}

    # FIMO columns: motif_id, motif_alt_id, sequence_name, start, stop, strand, score, p-value, q-value, matched_sequence
    # We need: sequence_name (anchor), start, stop, strand, p-value

    # Verify we have the required columns
    required_cols = ['sequence_name', 'start', 'stop', 'strand', 'p-value']
    missing_cols = [col for col in required_cols if col not in fimo_df.columns]
    if missing_cols:
        print(f"[motif_analysis][error] Missing columns in FIMO output: {missing_cols}")
        print(f"[motif_analysis][error] Available columns: {list(fimo_df.columns)}")
        return {}

    # Check for NaN values in required columns BEFORE processing
    nan_mask = fimo_df[required_cols].isnull().any(axis=1)
    if nan_mask.any():
        num_nan_rows = nan_mask.sum()
        print(f"[motif_analysis][warning] Found {num_nan_rows} rows with NaN values in required columns")
        print(f"[motif_analysis][warning] Removing these rows before parsing")
        fimo_df = fimo_df.dropna(subset=required_cols)

    if len(fimo_df) == 0:
        print(f"[motif_analysis][error] No valid FIMO results after filtering out NaN rows")
        return {}

    # Parse each FIMO hit
    try:
        for _, row in fimo_df.iterrows():
            anchor_id = row['sequence_name']
            motif_start = int(row['start'])
            motif_end = int(row['stop'])
            strand = row['strand']
            pvalue = float(row['p-value'])

            # Initialize anchor entry if not exists
            if anchor_id not in anchor_hits:
                anchor_hits[anchor_id] = {
                    'forward': None,
                    'reverse': None
                }

            # Store hit for this strand if it's better than existing
            strand_key = 'forward' if strand == '+' else 'reverse'

            # Only update if this is the best (lowest p-value) hit for this strand
            if anchor_hits[anchor_id][strand_key] is None or pvalue < anchor_hits[anchor_id][strand_key]['pvalue']:
                anchor_hits[anchor_id][strand_key] = {
                    'start': motif_start,
                    'end': motif_end,
                    'pvalue': pvalue,
                    'strand': strand
                }

    except (ValueError, TypeError) as e:
        print(f"[motif_analysis][error] Failed to parse FIMO row values: {e}")
        print(f"[motif_analysis][error] Check that start, stop, and p-value columns contain numeric values")
        return {}

    num_anchors_with_hits = sum(1 for hits in anchor_hits.values() if hits['forward'] is not None or hits['reverse'] is not None)
    print(f"[motif_analysis][info] Parsed FIMO results: {len(anchor_hits)} anchors, {num_anchors_with_hits} with hits")

    return anchor_hits


def determine_anchor_states(anchor_hits, bedpe_df, fdr_threshold=0.05, optimize_convergent=False):
    """
    Determine strand orientation state (F, R, or None) for each anchor.

    Two algorithms:
    1. Best q-value mode (default): Each anchor independently picks strand with lowest p-value
    2. Optimize convergent mode: Considers both anchors, prefers F-R convergent patterns

    Parameters:
    -----------
    anchor_hits : dict
        Parsed FIMO results from parse_fimo_results()
    bedpe_df : pd.DataFrame
        Original BEDPE DataFrame with loop information
    fdr_threshold : float, default=0.05
        FDR/p-value threshold for filtering hits
    optimize_convergent : bool, default=False
        If True, use convergent optimization algorithm; if False, use best q-value algorithm

    Returns:
    --------
    dict
        Dictionary mapping anchor IDs to their state:
        {
            'L1A1': 'F',  # Forward
            'L1A2': 'R',  # Reverse
            'L2A1': None,  # No hit
            ...
        }
    """
    anchor_states = {}

    if optimize_convergent:
        anchor_states = _determine_anchor_states_convergent(anchor_hits, bedpe_df, fdr_threshold)
    else:
        anchor_states = _determine_anchor_states_best_qvalue(anchor_hits, fdr_threshold)

    return anchor_states


def _determine_anchor_states_best_qvalue(anchor_hits, fdr_threshold=0.05):
    """
    Best q-value algorithm: Each anchor independently picks strand with lowest p-value.

    For each anchor:
    - If only forward hit passes FDR: state = 'F'
    - If only reverse hit passes FDR: state = 'R'
    - If both pass FDR: pick strand with lower p-value
    - If neither passes FDR: state = None

    Parameters:
    -----------
    anchor_hits : dict
        Parsed FIMO results
    fdr_threshold : float
        P-value threshold for filtering

    Returns:
    --------
    dict
        Dictionary mapping anchor IDs to their state (F, R, or None)
    """
    anchor_states = {}

    for anchor_id, hits in anchor_hits.items():
        forward_hit = hits['forward']
        reverse_hit = hits['reverse']

        # Filter by FDR threshold
        has_forward = forward_hit is not None and forward_hit['pvalue'] <= fdr_threshold
        has_reverse = reverse_hit is not None and reverse_hit['pvalue'] <= fdr_threshold

        if has_forward and has_reverse:
            # Both strands pass threshold - pick the one with lower p-value
            if forward_hit['pvalue'] <= reverse_hit['pvalue']:
                anchor_states[anchor_id] = 'F'
            else:
                anchor_states[anchor_id] = 'R'
        elif has_forward:
            anchor_states[anchor_id] = 'F'
        elif has_reverse:
            anchor_states[anchor_id] = 'R'
        else:
            anchor_states[anchor_id] = None

    return anchor_states


def _determine_anchor_states_convergent(anchor_hits, bedpe_df, fdr_threshold=0.05):
    """
    Convergent optimization algorithm: Considers both anchors, prefers F-R patterns.

    For each loop:
    1. If both anchors have hits on opposite strands (F-R), prefer that convergent pattern
    2. Otherwise, fall back to best q-value algorithm

    Parameters:
    -----------
    anchor_hits : dict
        Parsed FIMO results
    bedpe_df : pd.DataFrame
        Original BEDPE DataFrame with loop information
    fdr_threshold : float
        P-value threshold for filtering

    Returns:
    --------
    dict
        Dictionary mapping anchor IDs to their state (F, R, or None)
    """
    anchor_states = {}

    # First, get all hits that pass FDR threshold
    passing_hits = {}
    for anchor_id, hits in anchor_hits.items():
        passing_hits[anchor_id] = {}
        if hits['forward'] is not None and hits['forward']['pvalue'] <= fdr_threshold:
            passing_hits[anchor_id]['forward'] = hits['forward']
        if hits['reverse'] is not None and hits['reverse']['pvalue'] <= fdr_threshold:
            passing_hits[anchor_id]['reverse'] = hits['reverse']

    # Process each loop to check for convergent pattern
    for loop_id, (loop_idx, row) in enumerate(bedpe_df.iterrows(), start=1):
        # Use loop_id sequentially starting from 1, not the DataFrame index
        anchor1_id = f'L{loop_id}A1'
        anchor2_id = f'L{loop_id}A2'

        # Check if this loop has convergent F-R pattern
        has_convergent = (
            anchor1_id in passing_hits and 'forward' in passing_hits[anchor1_id] and
            anchor2_id in passing_hits and 'reverse' in passing_hits[anchor2_id]
        )

        if has_convergent:
            # Prefer convergent pattern
            anchor_states[anchor1_id] = 'F'
            anchor_states[anchor2_id] = 'R'
        else:
            # Fall back to best q-value for this anchor
            for anchor_id in [anchor1_id, anchor2_id]:
                if anchor_id not in anchor_states:  # Only process if not already assigned
                    if anchor_id in anchor_hits:
                        forward_hit = anchor_hits[anchor_id]['forward']
                        reverse_hit = anchor_hits[anchor_id]['reverse']

                        # Filter by FDR
                        has_forward = forward_hit is not None and forward_hit['pvalue'] <= fdr_threshold
                        has_reverse = reverse_hit is not None and reverse_hit['pvalue'] <= fdr_threshold

                        if has_forward and has_reverse:
                            if forward_hit['pvalue'] <= reverse_hit['pvalue']:
                                anchor_states[anchor_id] = 'F'
                            else:
                                anchor_states[anchor_id] = 'R'
                        elif has_forward:
                            anchor_states[anchor_id] = 'F'
                        elif has_reverse:
                            anchor_states[anchor_id] = 'R'
                        else:
                            anchor_states[anchor_id] = None

    return anchor_states


def classify_loops(bedpe_df, anchor_states):
    """
    Classify loops by anchor count and orientation pattern.

    Creates 7 categories:
    - zero_anchors: No CTCF binding at either anchor
    - one_anchor: CTCF binding at exactly one anchor
    - convergent: F-R orientation (both anchors, inward facing)
    - tandem_forward: F-F orientation (both anchors, same forward)
    - tandem_reverse: R-R orientation (both anchors, same reverse)
    - divergent: R-F orientation (both anchors, outward facing)
    - all_two_anchors: Combined category of all 2-anchor loops

    Parameters:
    -----------
    bedpe_df : pd.DataFrame
        Original BEDPE DataFrame with loop information
    anchor_states : dict
        Anchor states from determine_anchor_states()

    Returns:
    --------
    dict
        Dictionary with classifications:
        {
            'classifications': [category1, category2, ...],  # One per loop
            'category_counts': {
                'zero_anchors': int,
                'one_anchor': int,
                'convergent': int,
                'tandem_forward': int,
                'tandem_reverse': int,
                'divergent': int,
                'all_two_anchors': int
            }
        }
    """
    classifications = []
    category_counts = {
        'zero_anchors': 0,
        'one_anchor': 0,
        'convergent': 0,
        'tandem_forward': 0,
        'tandem_reverse': 0,
        'divergent': 0,
        'all_two_anchors': 0
    }

    for loop_id, (loop_idx, row) in enumerate(bedpe_df.iterrows(), start=1):
        # Use loop_id sequentially starting from 1, not the DataFrame index
        anchor1_id = f'L{loop_id}A1'
        anchor2_id = f'L{loop_id}A2'

        # Get states for both anchors (None if no CTCF binding)
        state1 = anchor_states.get(anchor1_id, None)
        state2 = anchor_states.get(anchor2_id, None)

        # Classify based on anchor states
        if state1 is None and state2 is None:
            # No CTCF at either anchor
            category = 'zero_anchors'
        elif state1 is None or state2 is None:
            # CTCF at exactly one anchor
            category = 'one_anchor'
        else:
            # CTCF at both anchors - determine orientation pattern
            orientation = f'{state1}-{state2}'

            if orientation == 'F-R':
                category = 'convergent'
            elif orientation == 'F-F':
                category = 'tandem_forward'
            elif orientation == 'R-R':
                category = 'tandem_reverse'
            elif orientation == 'R-F':
                category = 'divergent'
            else:
                # Should not happen, but handle gracefully
                category = 'all_two_anchors'

            # Also count in all_two_anchors
            category_counts['all_two_anchors'] += 1

        classifications.append(category)
        category_counts[category] += 1

    return {
        'classifications': classifications,
        'category_counts': category_counts
    }


def generate_statistics(bedpe_df, classifications, category_counts, output_dir):
    """
    Generate detailed orientation statistics and save to CSV.

    Calculates counts, percentages, and detailed breakdowns for all 7 categories.
    Saves statistics to CSV file and prints summary to console.

    Parameters:
    -----------
    bedpe_df : pd.DataFrame
        Original BEDPE DataFrame with loop information
    classifications : list
        Loop classifications from classify_loops()
    category_counts : dict
        Category counts from classify_loops()
    output_dir : str
        Output directory for statistics CSV file

    Returns:
    --------
    dict
        Statistics dictionary with all calculations
    """
    total_loops = len(bedpe_df)

    # Calculate percentages for each category
    stats = {
        'total_loops': total_loops,
        'zero_anchors': {
            'count': category_counts['zero_anchors'],
            'percentage': (category_counts['zero_anchors'] / total_loops * 100) if total_loops > 0 else 0
        },
        'one_anchor': {
            'count': category_counts['one_anchor'],
            'percentage': (category_counts['one_anchor'] / total_loops * 100) if total_loops > 0 else 0
        },
        'convergent': {
            'count': category_counts['convergent'],
            'percentage': (category_counts['convergent'] / total_loops * 100) if total_loops > 0 else 0
        },
        'tandem_forward': {
            'count': category_counts['tandem_forward'],
            'percentage': (category_counts['tandem_forward'] / total_loops * 100) if total_loops > 0 else 0
        },
        'tandem_reverse': {
            'count': category_counts['tandem_reverse'],
            'percentage': (category_counts['tandem_reverse'] / total_loops * 100) if total_loops > 0 else 0
        },
        'divergent': {
            'count': category_counts['divergent'],
            'percentage': (category_counts['divergent'] / total_loops * 100) if total_loops > 0 else 0
        },
        'all_two_anchors': {
            'count': category_counts['all_two_anchors'],
            'percentage': (category_counts['all_two_anchors'] / total_loops * 100) if total_loops > 0 else 0
        }
    }

    # Calculate percentages within 2-anchor loops only
    two_anchor_total = category_counts['convergent'] + category_counts['tandem_forward'] + category_counts['tandem_reverse'] + category_counts['divergent']

    stats['convergent']['percentage_of_2anchor'] = (category_counts['convergent'] / two_anchor_total * 100) if two_anchor_total > 0 else 0
    stats['tandem_forward']['percentage_of_2anchor'] = (category_counts['tandem_forward'] / two_anchor_total * 100) if two_anchor_total > 0 else 0
    stats['tandem_reverse']['percentage_of_2anchor'] = (category_counts['tandem_reverse'] / two_anchor_total * 100) if two_anchor_total > 0 else 0
    stats['divergent']['percentage_of_2anchor'] = (category_counts['divergent'] / two_anchor_total * 100) if two_anchor_total > 0 else 0

    # Print summary to console
    print("\n[motif_analysis][info] ===== Loop Classification Statistics =====")
    print(f"[motif_analysis][info] Total loops analyzed: {total_loops}")
    print(f"[motif_analysis][info]")
    print(f"[motif_analysis][info] By anchor count:")
    print(f"[motif_analysis][info]   Zero anchors:  {stats['zero_anchors']['count']:>8} ({stats['zero_anchors']['percentage']:>6.2f}%)")
    print(f"[motif_analysis][info]   One anchor:    {stats['one_anchor']['count']:>8} ({stats['one_anchor']['percentage']:>6.2f}%)")
    print(f"[motif_analysis][info]   Two anchors:   {stats['all_two_anchors']['count']:>8} ({stats['all_two_anchors']['percentage']:>6.2f}%)")
    print(f"[motif_analysis][info]")
    print(f"[motif_analysis][info] Two-anchor loop orientation breakdown:")
    print(f"[motif_analysis][info]   Convergent (F-R):      {stats['convergent']['count']:>8} ({stats['convergent']['percentage_of_2anchor']:>6.2f}% of 2-anchor)")
    print(f"[motif_analysis][info]   Tandem Forward (F-F):  {stats['tandem_forward']['count']:>8} ({stats['tandem_forward']['percentage_of_2anchor']:>6.2f}% of 2-anchor)")
    print(f"[motif_analysis][info]   Tandem Reverse (R-R):  {stats['tandem_reverse']['count']:>8} ({stats['tandem_reverse']['percentage_of_2anchor']:>6.2f}% of 2-anchor)")
    print(f"[motif_analysis][info]   Divergent (R-F):       {stats['divergent']['count']:>8} ({stats['divergent']['percentage_of_2anchor']:>6.2f}% of 2-anchor)")
    print(f"[motif_analysis][info] ==========================================")

    # Save to CSV
    os.makedirs(output_dir, exist_ok=True)
    csv_path = os.path.join(output_dir, 'loop_classification_statistics.csv')

    try:
        stats_data = {
            'Category': [
                'Zero_Anchors',
                'One_Anchor',
                'Convergent_FR',
                'Tandem_Forward_FF',
                'Tandem_Reverse_RR',
                'Divergent_RF',
                'All_Two_Anchors'
            ],
            'Count': [
                stats['zero_anchors']['count'],
                stats['one_anchor']['count'],
                stats['convergent']['count'],
                stats['tandem_forward']['count'],
                stats['tandem_reverse']['count'],
                stats['divergent']['count'],
                stats['all_two_anchors']['count']
            ],
            'Percentage_of_Total': [
                f"{stats['zero_anchors']['percentage']:.2f}",
                f"{stats['one_anchor']['percentage']:.2f}",
                f"{stats['convergent']['percentage']:.2f}",
                f"{stats['tandem_forward']['percentage']:.2f}",
                f"{stats['tandem_reverse']['percentage']:.2f}",
                f"{stats['divergent']['percentage']:.2f}",
                f"{stats['all_two_anchors']['percentage']:.2f}"
            ],
            'Percentage_of_2Anchor': [
                'N/A',
                'N/A',
                f"{stats['convergent']['percentage_of_2anchor']:.2f}",
                f"{stats['tandem_forward']['percentage_of_2anchor']:.2f}",
                f"{stats['tandem_reverse']['percentage_of_2anchor']:.2f}",
                f"{stats['divergent']['percentage_of_2anchor']:.2f}",
                'N/A'
            ]
        }

        stats_df = pd.DataFrame(stats_data)
        stats_df.to_csv(csv_path, index=False)
        print(f"[motif_analysis][info] Statistics saved to: {csv_path}")

    except Exception as e:
        print(f"[motif_analysis][error] Failed to save statistics CSV: {e}")

    return stats


def count_hits_per_anchor(anchor_all_hits, bedpe_df):
    """
    Count the number of CTCF hits per anchor across all loops.

    Parameters:
    -----------
    anchor_all_hits : dict
        Dictionary mapping anchor IDs to list of all hits (from parse_fimo_results_all_hits)
    bedpe_df : pd.DataFrame
        Original BEDPE DataFrame with loop information

    Returns:
    --------
    dict
        Dictionary with hit count statistics:
        {
            'hit_counts': [0, 1, 2, 3, ...],  # List of hit counts per anchor
            'count_distribution': {0: n0, 1: n1, 2: n2, ...},  # Frequency of each count
            'total_anchors': int,
            'anchors_with_hits': int,
            'mean_hits': float,
            'median_hits': float
        }
    """
    import numpy as np

    # Calculate total number of anchors (2 per loop)
    total_anchors = len(bedpe_df) * 2

    # Count hits per anchor
    hit_counts = []
    for loop_id in range(1, len(bedpe_df) + 1):
        anchor1_id = f'L{loop_id}A1'
        anchor2_id = f'L{loop_id}A2'

        # Count hits for each anchor
        count1 = len(anchor_all_hits.get(anchor1_id, []))
        count2 = len(anchor_all_hits.get(anchor2_id, []))

        hit_counts.append(count1)
        hit_counts.append(count2)

    # Calculate distribution
    count_distribution = {}
    for count in hit_counts:
        count_distribution[count] = count_distribution.get(count, 0) + 1

    # Calculate statistics
    anchors_with_hits = sum(1 for c in hit_counts if c > 0)
    mean_hits = np.mean(hit_counts)
    median_hits = np.median(hit_counts)

    # Print summary
    print(f"\n[motif_analysis][info] ===== CTCF Hit Distribution =====")
    print(f"[motif_analysis][info] Total anchors: {total_anchors}")
    print(f"[motif_analysis][info] Anchors with hits: {anchors_with_hits} ({anchors_with_hits/total_anchors*100:.1f}%)")
    print(f"[motif_analysis][info] Mean hits per anchor: {mean_hits:.2f}")
    print(f"[motif_analysis][info] Median hits per anchor: {median_hits:.1f}")
    print(f"[motif_analysis][info]")
    print(f"[motif_analysis][info] Hit count distribution:")

    # Sort by hit count for display
    for count in sorted(count_distribution.keys()):
        num_anchors = count_distribution[count]
        pct = num_anchors / total_anchors * 100
        print(f"[motif_analysis][info]   {count} hits: {num_anchors:>8} anchors ({pct:>5.1f}%)")

    print(f"[motif_analysis][info] =====================================\n")

    return {
        'hit_counts': hit_counts,
        'count_distribution': count_distribution,
        'total_anchors': total_anchors,
        'anchors_with_hits': anchors_with_hits,
        'mean_hits': mean_hits,
        'median_hits': median_hits
    }


def save_hit_distribution_stats(hit_stats, output_dir):
    """
    Save CTCF hit distribution statistics to CSV.

    Parameters:
    -----------
    hit_stats : dict
        Hit count statistics from count_hits_per_anchor()
    output_dir : str
        Output directory for CSV file

    Returns:
    --------
    str
        Path to saved CSV file
    """
    import os

    os.makedirs(output_dir, exist_ok=True)

    # Create DataFrame from count distribution
    count_dist = hit_stats['count_distribution']
    total_anchors = hit_stats['total_anchors']

    data = []
    for count in sorted(count_dist.keys()):
        num_anchors = count_dist[count]
        pct = num_anchors / total_anchors * 100
        data.append({
            'num_hits': count,
            'num_anchors': num_anchors,
            'percentage': f"{pct:.2f}"
        })

    df = pd.DataFrame(data)
    csv_path = os.path.join(output_dir, 'ctcf_hit_distribution.csv')
    df.to_csv(csv_path, index=False)

    print(f"[motif_analysis][info] CTCF hit distribution saved to: {csv_path}")

    return csv_path


def output_separated_bedpe(bedpe_df, classifications, anchor_states, output_dir):
    """
    Output separated BEDPE files for each classification category and a comprehensive file.

    Creates 7 separate BEDPE files:
    - zero_anchors.bedpe
    - one_anchor.bedpe
    - convergent.bedpe
    - tandem_forward.bedpe
    - tandem_reverse.bedpe
    - divergent.bedpe
    - all_two_anchors.bedpe

    Also creates a comprehensive TSV file with all loops and their classifications:
    - all_loops_with_classifications.tsv

    Each file includes original BEDPE columns plus anchor state information.

    Parameters:
    -----------
    bedpe_df : pd.DataFrame
        Original BEDPE DataFrame with loop information
    classifications : list
        Loop classifications from classify_loops()
    anchor_states : dict
        Anchor states from determine_anchor_states()
    output_dir : str
        Output directory for output files

    Returns:
    --------
    dict
        Dictionary with file paths for each category and comprehensive file
    """
    os.makedirs(output_dir, exist_ok=True)

    # Initialize output dataframes for each category
    output_dfs = {
        'zero_anchors': [],
        'one_anchor': [],
        'convergent': [],
        'tandem_forward': [],
        'tandem_reverse': [],
        'divergent': [],
        'all_two_anchors': []
    }

    # Data for comprehensive file
    comprehensive_data = []

    # Process each loop - single pass
    for loop_pos, (loop_idx, row) in enumerate(bedpe_df.iterrows()):
        loop_id = loop_pos + 1
        # Use loop_pos sequentially (0-based) to index classifications list
        category = classifications[loop_pos]
        anchor1_id = f'L{loop_id}A1'
        anchor2_id = f'L{loop_id}A2'

        # Get anchor states
        state1 = anchor_states.get(anchor1_id, None)
        state2 = anchor_states.get(anchor2_id, None)

        # Create output row with anchor state information
        output_row = row.copy()
        output_row['anchor1_state'] = state1 if state1 else 'None'
        output_row['anchor2_state'] = state2 if state2 else 'None'

        # Add to appropriate category
        output_dfs[category].append(output_row)

        # Also add 2-anchor loops to all_two_anchors
        if category in ['convergent', 'tandem_forward', 'tandem_reverse', 'divergent']:
            output_dfs['all_two_anchors'].append(output_row)

        # Build record for comprehensive file
        orientation = f'{state1}-{state2}' if state1 and state2 else 'N/A'
        comprehensive_record = {
            'loop_id': loop_id,
            'chr1': row['chr1'],
            'start1': row['start1'],
            'end1': row['end1'],
            'chr2': row['chr2'],
            'start2': row['start2'],
            'end2': row['end2'],
            'anchor1_state': state1 if state1 else 'None',
            'anchor2_state': state2 if state2 else 'None',
            'orientation': orientation,
            'classification': category
        }

        # Add any additional columns from original BEDPE
        for col in row.index:
            if col not in ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']:
                comprehensive_record[col] = row[col]

        comprehensive_data.append(comprehensive_record)

    # Write each category to file
    file_paths = {}
    categories = ['zero_anchors', 'one_anchor', 'convergent', 'tandem_forward', 'tandem_reverse', 'divergent', 'all_two_anchors']

    for category in categories:
        if len(output_dfs[category]) > 0:
            # Convert list of rows back to DataFrame
            category_df = pd.DataFrame(output_dfs[category])

            # Write to file
            file_path = os.path.join(output_dir, f'{category}.bedpe')
            category_df.to_csv(file_path, sep='\t', index=False)
            file_paths[category] = file_path
            print(f"[motif_analysis][info] Saved {category}: {len(category_df)} loops to {file_path}")
        else:
            file_path = os.path.join(output_dir, f'{category}.bedpe')
            file_paths[category] = file_path
            print(f"[motif_analysis][info] No loops in category '{category}' - empty file created")

    # Write comprehensive loops file
    comprehensive_df = pd.DataFrame(comprehensive_data)
    column_order = [
        'loop_id', 'chr1', 'start1', 'end1', 'chr2', 'start2', 'end2',
        'anchor1_state', 'anchor2_state', 'orientation', 'classification'
    ]
    # Add any remaining columns
    remaining_cols = [col for col in comprehensive_df.columns if col not in column_order]
    column_order.extend(remaining_cols)
    comprehensive_df = comprehensive_df[column_order]

    comprehensive_path = os.path.join(output_dir, 'all_loops_with_classifications.tsv')
    comprehensive_df.to_csv(comprehensive_path, sep='\t', index=False)
    file_paths['comprehensive'] = comprehensive_path
    print(f"[motif_analysis][info] Saved comprehensive loops file: {comprehensive_path}")

    return file_paths
