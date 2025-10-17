"""
APA (Aggregate Peak Analysis) using cooltools pileup.

Aggregates Hi-C contact frequencies around loop anchors to assess loop strength.
"""

import os
import subprocess
import numpy as np
import pandas as pd
import cooler
from cooltools import pileup
import warnings

from .utils import calculate_apa_statistics
from .visualizations import plot_apa_comparison_heatmaps

# Suppress pandas FutureWarnings from cooltools (until cooltools updates)
warnings.filterwarnings('ignore', category=FutureWarning, module='cooltools')


def balance_cooler_matrix(hic_file, nproc=1):
    """
    Balance a cooler Hi-C matrix using iterative correction.

    Parameters:
    -----------
    hic_file : str
        Path to Hi-C cooler file
    nproc : int
        Number of threads to use

    Returns:
    --------
    bool
        True if balancing succeeded, False otherwise
    """
    try:
        print(f"[apa][info] Balancing Hi-C matrix: {hic_file}")
        result = subprocess.run(
            ['cooler', 'balance', hic_file, '--nproc', str(nproc)],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            print(f"[apa][info] Matrix balanced successfully")
            return True
        else:
            print(f"[apa][error] Cooler balance failed: {result.stderr}")
            return False
    except Exception as e:
        print(f"[apa][error] Failed to balance matrix: {e}")
        return False


def ensure_balanced_matrix(hic_file, resolution=None, nproc=1):
    """
    Check if Hi-C matrix is balanced. If not, attempt to balance it.

    Parameters:
    -----------
    hic_file : str
        Path to Hi-C cooler file (can be mcool with multiple resolutions)
    resolution : int, optional
        Resolution in bp (for mcool files with multiple resolutions)
    nproc : int
        Number of threads to use for balancing

    Returns:
    --------
    tuple
        (hic_cooler_object, is_balanced_bool)
    """
    # If it's an mcool file and resolution is specified, construct the group path
    hic_file_to_load = hic_file
    if resolution and hic_file.endswith('.cool') and '::' not in hic_file:
        # For mcool files, append the resolution group path
        hic_file_to_load = f"{hic_file}::/resolutions/{resolution}"

    try:
        hic = cooler.Cooler(hic_file_to_load)
        if '::' in hic_file_to_load:
            print(f"[apa][info] Loaded mcool file at {resolution}bp resolution")
    except Exception as e:
        print(f"[apa][error] Failed to load Hi-C file: {e}")
        return None, False

    # Check if already balanced
    if 'weight' in hic.bins().columns:
        print(f"[apa][info] Hi-C matrix is already balanced")
        return hic, True

    print(f"[apa][warning] Hi-C matrix is not balanced (missing 'weight' column)")
    print(f"[apa][info] Attempting to balance matrix...")

    # Try to balance the matrix (balance the base file, not the resolution group)
    if balance_cooler_matrix(hic_file, nproc=nproc):
        # Reload the cooler file to get the balanced version
        try:
            hic = cooler.Cooler(hic_file_to_load)
            if 'weight' in hic.bins().columns:
                print(f"[apa][info] Matrix is now balanced and ready for APA")
                return hic, True
        except Exception as e:
            print(f"[apa][error] Failed to reload balanced matrix: {e}")

    print(f"[apa][error] Could not balance Hi-C matrix")
    return None, False


def run_apa_analysis(loops_df, apa_config, output_dir):
    """
    Run APA analysis on loops using cooltools pileup.

    Parameters:
    -----------
    loops_df : pd.DataFrame
        BEDPE dataframe with columns: chr1, start1, end1, chr2, start2, end2
    apa_config : dict
        APA configuration with keys: hic_file, window_size, resolution
    output_dir : str
        Output directory for results

    Returns:
    --------
    dict
        APA results dictionary with aggregated matrices and statistics
    """
    hic_file = apa_config['hic_file']
    window_size = apa_config['window_size']
    resolution = apa_config['resolution']
    nproc = apa_config.get('nproc', 1)

    print(f"[apa][info] Running APA analysis on {len(loops_df)} loops")
    print(f"[apa][info] Hi-C file: {hic_file}")
    print(f"[apa][info] Window size: ±{window_size/1000:.0f}kb, Resolution: {resolution/1000:.0f}kb")
    print(f"[apa][info] Using {nproc} thread(s)")

    # Check and balance Hi-C matrix if needed
    hic, is_balanced = ensure_balanced_matrix(hic_file, resolution=resolution, nproc=nproc)

    if not is_balanced or hic is None:
        print(f"[apa][error] Cannot proceed without a balanced Hi-C matrix")
        return None

    try:
        print(f"[apa][info] Loaded Hi-C matrix: {hic.info}")
    except Exception as e:
        print(f"[apa][error] Failed to access Hi-C matrix: {e}")
        return None

    # Convert loop coordinates to regions for pileup
    # pileup expects a DataFrame with columns: chrom1, start1, end1, chrom2, start2, end2
    regions_data = []
    for idx, row in loops_df.iterrows():
        # Normalize chromosome names to match Hi-C file format
        chr1 = row['chr1']
        chr2 = row['chr2']

        # Add 'chr' prefix if not present
        if not chr1.startswith('chr'):
            chr1 = f'chr{chr1}'
        if not chr2.startswith('chr'):
            chr2 = f'chr{chr2}'

        # Region around first anchor
        anchor1_mid = (row['start1'] + row['end1']) // 2
        region1_start = max(0, anchor1_mid - window_size // 2)
        region1_end = anchor1_mid + window_size // 2

        # Region around second anchor
        anchor2_mid = (row['start2'] + row['end2']) // 2
        region2_start = max(0, anchor2_mid - window_size // 2)
        region2_end = anchor2_mid + window_size // 2

        regions_data.append({
            'chrom1': chr1,
            'start1': region1_start,
            'end1': region1_end,
            'chrom2': chr2,
            'start2': region2_start,
            'end2': region2_end
        })

    regions_df = pd.DataFrame(regions_data)

    # Run pileup
    try:
        pileup_results = pileup(hic, regions_df, nproc=nproc)
        print(f"[apa][info] Pileup completed successfully")
    except Exception as e:
        print(f"[apa][error] Pileup failed: {e}")
        return None

    # Aggregate across all loops (sum along the first dimension)
    # IMPORTANT: The pileup results are asymmetric due to the Hi-C matrix storage format (symmetric-upper).
    #
    # Root cause: Hi-C matrices in balanced cooler files use 'symmetric-upper' storage, where:
    #   - Lower triangle: contains the actual contact frequencies
    #   - Upper triangle: is a mirror of the lower triangle
    #
    # When cooltools.pileup() extracts regions, it treats region1 as rows and region2 as columns.
    # For intra-chromosomal loops where pos1 < pos2, the extracted matrix is always from the
    # upper triangle (which has less signal than the lower triangle), causing systematic bias.
    #
    # Solution: Symmetrize the aggregated matrix by averaging it with its transpose.
    # This effectively treats the loop in both orientations (pos1→pos2 and pos2→pos1),
    # resulting in a balanced APA visualization.
    aggregated = np.nansum(pileup_results, axis=0)
    aggregated_symmetric = (aggregated + aggregated.T) / 2
    print(f"[apa][info] Symmetrized aggregated matrix to correct for storage format bias")

    # Generate APA statistics using the symmetrized matrix
    apa_stats = calculate_apa_statistics(aggregated_symmetric)

    # Save aggregated matrix (use symmetrized version)
    np.save(os.path.join(output_dir, 'apa_aggregated_matrix.npy'), aggregated_symmetric)

    # Save statistics
    stats_df = pd.DataFrame([apa_stats])
    stats_df.to_csv(os.path.join(output_dir, 'apa_statistics.csv'), index=False)

    print(f"[apa][info] APA results saved to {output_dir}")

    return {
        'aggregated_matrix': aggregated_symmetric,
        'statistics': apa_stats,
        'output_dir': output_dir
    }


def run_apa_by_category(categories_dict, apa_config, output_dir):
    """
    Run APA analysis on loops grouped by category.

    Parameters:
    -----------
    categories_dict : dict
        Dictionary mapping category names to loop DataFrames
    apa_config : dict
        APA configuration
    output_dir : str
        Output directory for results

    Returns:
    --------
    dict
        APA results by category
    """
    results = {}

    for category, loops_df in categories_dict.items():
        if len(loops_df) == 0:
            print(f"[apa][warning] Skipping empty category: {category}")
            continue

        print(f"[apa][info] Running APA for category: {category} ({len(loops_df)} loops)")

        category_output_dir = os.path.join(output_dir, f'apa_results/{category}')
        os.makedirs(category_output_dir, exist_ok=True)

        result = run_apa_analysis(loops_df, apa_config, category_output_dir)
        results[category] = result

    return results
