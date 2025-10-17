import pandas as pd
import numpy as np
import os


def read_bedpe(filepath, has_header=True):
    """
    Read a BEDPE file and return as a pandas DataFrame.

    Parameters:
    -----------
    filepath : str
        Path to the BEDPE file
    has_header : bool, default=True
        Whether the BEDPE file has a header row

    Returns:
    --------
    pd.DataFrame
        DataFrame containing the BEDPE data with columns:
        chr1, start1, end1, chr2, start2, end2, [optional additional columns]
    """
    # Define expected column names for BEDPE format
    column_names = ['chr1', 'start1', 'end1', 'chr2', 'start2', 'end2']

    # Read the BEDPE file
    if has_header:
        # Read with header and ensure we have the expected columns
        df = pd.read_csv(filepath, sep='\t')
        # Rename columns to standardized names if they don't already match
        if len(df.columns) >= 6:
            # Check if columns already match standard names
            if list(df.columns[:6]) != column_names:
                # Rename the first 6 columns to standard names
                df.columns = column_names + list(df.columns[6:])
                print(f"[utils][info] Renamed columns to standard BEDPE format")
    else:
        # Read without header
        df = pd.read_csv(filepath, sep='\t', header=None)
        # Assign standard column names
        num_cols = len(df.columns)
        df.columns = column_names + [f'col_{i}' for i in range(6, num_cols)]

    # Print number of loops
    num_loops = len(df)
    print(f"[utils][info] Successfully read BEDPE file: {filepath}")
    print(f"[utils][info] Number of loops: {num_loops}")

    return df


def filter_intrachromosomal(df):
    """
    Filter loops to only include intra-chromosomal interactions (same chromosome).

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with columns: chr1, chr2

    Returns:
    --------
    pd.DataFrame
        Filtered DataFrame containing only intra-chromosomal loops
    """
    intra_df = df[df['chr1'] == df['chr2']].copy()
    num_filtered = len(df) - len(intra_df)
    print(f"[utils][info] Filtered to intra-chromosomal loops: {len(intra_df)} loops")
    print(f"[utils][info] Removed {num_filtered} inter-chromosomal loops")

    return intra_df


def filter_interchromosomal(df):
    """
    Filter loops to only include inter-chromosomal interactions (different chromosomes).

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with columns: chr1, chr2

    Returns:
    --------
    pd.DataFrame
        Filtered DataFrame containing only inter-chromosomal loops
    """
    inter_df = df[df['chr1'] != df['chr2']].copy()
    num_filtered = len(df) - len(inter_df)
    print(f"[utils][info] Filtered to inter-chromosomal loops: {len(inter_df)} loops")
    print(f"[utils][info] Removed {num_filtered} intra-chromosomal loops")

    return inter_df


def filter_by_loop_length(df, min_length_bp):
    """
    Filter loops to only include those with length >= min_length_bp.

    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame with columns: start1, end1, start2, end2
    min_length_bp : int
        Minimum loop length in base pairs

    Returns:
    --------
    pd.DataFrame
        Filtered DataFrame containing only loops meeting minimum length criteria
    """
    # Calculate loop lengths (distance between anchor midpoints)
    anchor1_mid = (df['start1'] + df['end1']) / 2
    anchor2_mid = (df['start2'] + df['end2']) / 2
    loop_lengths = np.abs(anchor2_mid - anchor1_mid)

    # Filter
    filtered_df = df[loop_lengths >= min_length_bp].copy()
    num_removed = len(df) - len(filtered_df)
    pct_kept = (len(filtered_df) / len(df) * 100) if len(df) > 0 else 0

    return filtered_df


def setup_output_directories(output_base):
    """
    Create output directory structure with subdirectories for temp, figures, and generated files.

    Parameters:
    -----------
    output_base : str
        Base path for output directory

    Returns:
    --------
    dict
        Dictionary containing paths to subdirectories:
        {'temp': path, 'figures': path, 'generated_files': path}
    """
    # Create main output directory
    os.makedirs(output_base, exist_ok=True)

    # Create subdirectories
    temp_dir = os.path.join(output_base, 'temp')
    figures_dir = os.path.join(output_base, 'figures')
    generated_files_dir = os.path.join(output_base, 'generated_files')

    os.makedirs(temp_dir, exist_ok=True)
    os.makedirs(figures_dir, exist_ok=True)
    os.makedirs(generated_files_dir, exist_ok=True)

    print(f"[utils][info] Output directory structure created at: {output_base}")
    print(f"[utils][info]   - temp: {temp_dir}")
    print(f"[utils][info]   - figures: {figures_dir}")
    print(f"[utils][info]   - generated_files: {generated_files_dir}")

    return {
        'temp': temp_dir,
        'figures': figures_dir,
        'generated_files': generated_files_dir,
        'base': output_base
    }


def calculate_apa_statistics(aggregated_matrix):
    """
    Calculate statistics from aggregated APA matrix.

    Parameters:
    -----------
    aggregated_matrix : np.ndarray
        Aggregated contact matrix (2D)

    Returns:
    --------
    dict
        Statistics including corner values and enrichment
    """
    if aggregated_matrix is None or aggregated_matrix.size == 0:
        return {}

    matrix = aggregated_matrix

    # Check if matrix is valid 2D
    if matrix.ndim != 2:
        print(f"[utils][warning] Expected 2D matrix, got shape: {matrix.shape}")
        return {}

    if matrix.shape[0] < 2 or matrix.shape[1] < 2:
        print(f"[utils][warning] Matrix too small for statistics: {matrix.shape}")
        return {}

    center_idx_x = matrix.shape[0] // 2
    center_idx_y = matrix.shape[1] // 2

    # Extract corner values
    top_left = np.nanmean(matrix[:center_idx_x, :center_idx_y])
    top_right = np.nanmean(matrix[:center_idx_x, center_idx_y:])
    bottom_left = np.nanmean(matrix[center_idx_x:, :center_idx_y])
    bottom_right = np.nanmean(matrix[center_idx_x:, center_idx_y:])

    # Center value (peak)
    center_value = matrix[center_idx_x, center_idx_y]

    # Calculate enrichment relative to background (corners)
    background = np.nanmean([top_left, top_right, bottom_left, bottom_right])
    enrichment = center_value / background if background > 0 else 0

    return {
        'center_value': float(center_value) if not np.isnan(center_value) else 0,
        'background_mean': float(background) if not np.isnan(background) else 0,
        'enrichment': float(enrichment) if not np.isnan(enrichment) else 0,
        'top_left': float(top_left) if not np.isnan(top_left) else 0,
        'top_right': float(top_right) if not np.isnan(top_right) else 0,
        'bottom_left': float(bottom_left) if not np.isnan(bottom_left) else 0,
        'bottom_right': float(bottom_right) if not np.isnan(bottom_right) else 0
    }

