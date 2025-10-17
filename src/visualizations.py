import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from scipy import stats

# Hard-coded figure parameters
FONT_SIZE = 18
VIOLIN_FIGURE_WIDTH = 16  # Wider for violin plots with multiple categories
VIOLIN_FIGURE_HEIGHT = 6
APA_FIGURE_WIDTH = 8  # Per subplot
APA_FIGURE_HEIGHT = 7  # Per subplot


def plot_loop_length_violin_by_category(bedpe_dfs, figure_dir, generated_files_dir=None, filename='loop_length_violin_by_category.png'):
    """
    Create a violin plot showing loop length distributions across all categories.

    Displays distributions in order:
    1. All loops
    2. Zero motifs
    3. One motif
    4. All two motifs
    5. Convergent (F-R)
    6. Tandem forward (F-F)
    7. Tandem reverse (R-R)
    8. Divergent (R-F)

    Parameters:
    -----------
    bedpe_dfs : dict
        Dictionary mapping category names to DataFrames with loop information.
        Expected keys: 'all_loops', 'zero_anchors', 'one_anchor', 'all_two_anchors',
                      'convergent', 'tandem_forward', 'tandem_reverse', 'divergent'
    figure_dir : str
        Directory to save the output figure
    generated_files_dir : str, optional
        Directory to save statistics CSV. If None, saves to figure_dir
    filename : str, default='loop_length_violin_by_category.png'
        Name of the output file

    Returns:
    --------
    None
        Saves figure to figure_dir/filename and statistics to generated_files_dir/loop_length_statistics.csv
    """
    # Use figure_dir for stats if generated_files_dir not provided
    if generated_files_dir is None:
        generated_files_dir = figure_dir
    # Define category order and labels
    category_order = [
        'all_loops',
        'zero_anchors',
        'one_anchor',
        'all_two_anchors',
        'convergent',
        'tandem_forward',
        'tandem_reverse',
        'divergent'
    ]

    # Prepare data for violin plot
    plot_data = []
    category_labels = []
    category_data = {}  # Store category data for statistical testing

    for category in category_order:
        if category not in bedpe_dfs or len(bedpe_dfs[category]) == 0:
            continue

        df = bedpe_dfs[category]

        # Calculate loop lengths if not already present
        if 'loop_length' not in df.columns:
            # Calculate from anchor positions
            loop_lengths = np.abs(
                ((df['start2'] + df['end2']) / 2) - ((df['start1'] + df['end1']) / 2)
            )
        else:
            loop_lengths = df['loop_length'].values

        # Filter out zero lengths
        valid_lengths = loop_lengths[loop_lengths > 0]

        if len(valid_lengths) == 0:
            continue

        # Convert to Kbp
        lengths_kbp = valid_lengths / 1000

        # Store for statistical testing
        category_data[category] = lengths_kbp.values

        # Add to plot data
        for length in lengths_kbp:
            plot_data.append({
                'loop_length': length,
                'category': category
            })

        category_labels.append(category)

    if len(plot_data) == 0:
        print("[visualizations][warning] No valid loop lengths to plot")
        return

    # Convert to DataFrame for seaborn-style plotting
    plot_df = pd.DataFrame(plot_data)

    # Perform pairwise t-tests to identify significant differences
    # Compare each category to all_loops (reference)
    significance_markers = {}
    stats_data = []  # For CSV output

    if 'all_loops' in category_data:
        reference_data = category_data['all_loops']
        median_ref = np.median(reference_data)

        for cat_idx, category in enumerate(category_labels):
            if category == 'all_loops':
                continue

            cat_data = category_data[category]
            # Perform independent samples t-test
            t_stat, p_value = stats.ttest_ind(reference_data, cat_data)

            # Calculate fold change
            median_cat = np.median(cat_data)
            fold_change = median_cat / median_ref if median_ref > 0 else 0

            # Determine stars based on p-value
            # Standard significance levels: *** p<0.001, ** p<0.01, * p<0.05
            if p_value < 0.001:
                stars = '***'
            elif p_value < 0.01:
                stars = '**'
            elif p_value < 0.05:
                stars = '*'
            else:
                stars = ''

            # Store marker for plot (only if significant at p<0.05)
            if p_value < 0.05:
                significance_markers[cat_idx] = stars

            # Store statistics for CSV
            stats_data.append({
                'category': category,
                'p_value': p_value,
                'fold_change': fold_change,
                'median_ref': median_ref,
                'median_cat': median_cat,
                'significance': stars
            })

            print(f"[visualizations][info] {category} vs all_loops: p={p_value:.4e}, fold_change={fold_change:.2f}, significance={stars}")

    # Create figure and axis with wider dimensions
    fig, ax = plt.subplots(figsize=(VIOLIN_FIGURE_WIDTH, VIOLIN_FIGURE_HEIGHT))

    # Create violin plot
    parts = ax.violinplot(
        [plot_df[plot_df['category'] == cat]['loop_length'].values for cat in category_labels],
        positions=range(len(category_labels)),
        widths=0.7,
        showmeans=True,
        showmedians=True
    )

    # Customize violin plot colors
    for pc in parts['bodies']:
        pc.set_facecolor('steelblue')
        pc.set_alpha(0.7)
        pc.set_edgecolor('black')

    # Customize other elements
    for partname in ('cbars', 'cmins', 'cmaxes', 'cmedians', 'cmeans'):
        if partname in parts:
            parts[partname].set_edgecolor('black')
            parts[partname].set_linewidth(1.5)

    # Set x-axis labels and positions with significance markers
    ax.set_xticks(range(len(category_labels)))

    # Add significance markers to labels
    labeled_with_markers = []
    for idx, label in enumerate(category_labels):
        if idx in significance_markers:
            labeled_with_markers.append(f"{label}\n{significance_markers[idx]}")
        else:
            labeled_with_markers.append(label)

    ax.set_xticklabels(labeled_with_markers, rotation=45, ha='right', fontsize=FONT_SIZE)
    ax.tick_params(axis='y', labelsize=FONT_SIZE)

    # Labels
    ax.set_ylabel('Loop Length (Kbp)', fontsize=FONT_SIZE)

    # Grid for readability
    ax.yaxis.grid(True, alpha=0.3)
    ax.set_axisbelow(True)

    # Tight layout
    plt.tight_layout()

    # Save figure
    output_path = os.path.join(figure_dir, filename)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"[visualizations][info] Saved violin plot: {output_path}")

    # Save statistics to CSV
    if stats_data:
        stats_df = pd.DataFrame(stats_data)
        stats_csv_path = os.path.join(generated_files_dir, 'loop_length_statistics.csv')
        stats_df.to_csv(stats_csv_path, index=False)
        print(f"[visualizations][info] Saved loop length statistics: {stats_csv_path}")

    # Close the figure
    plt.close()


def plot_apa_comparison_heatmaps(categories_results, output_dir, resolution, window_size):
    """
    Create comparison APA heatmaps for all loop categories in one figure.

    Parameters:
    -----------
    categories_results : dict
        Dictionary mapping category names to APA results (with 'aggregated_matrix' key)
    output_dir : str
        Output directory for figure
    resolution : int
        Hi-C resolution in bp
    window_size : int
        Window size in bp
    """
    # Filter out None/empty results
    valid_results = {k: v for k, v in categories_results.items() if v and v.get('aggregated_matrix') is not None}

    if not valid_results:
        print(f"[visualizations][warning] No valid APA results to visualize")
        return

    n_categories = len(valid_results)
    n_cols = min(3, n_categories)
    n_rows = (n_categories + n_cols - 1) // n_cols

    # Increase figure size with extra space for colorbar
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(APA_FIGURE_WIDTH*n_cols, APA_FIGURE_HEIGHT*n_rows))
    axes = axes.flatten() if n_categories > 1 else [axes]

    # Determine global scale for consistent coloring
    all_matrices = [v['aggregated_matrix'] for v in valid_results.values()]
    vmax = np.nanpercentile([m.max() for m in all_matrices], 95)

    for idx, (category, result) in enumerate(valid_results.items()):
        ax = axes[idx]
        matrix = result['aggregated_matrix']
        matrix_log = np.log10(matrix + 1)

        im = ax.imshow(matrix_log, cmap='YlOrRd', aspect='auto', origin='lower', vmax=np.log10(vmax + 1))

        # Set axis labels (in kb)
        n_bins = matrix.shape[0]
        tick_positions = np.linspace(0, n_bins - 1, 5)
        tick_labels = [f"{-window_size//2 + int(i*window_size/(n_bins-1))/1000:.0f}" for i in range(5)]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=FONT_SIZE)
        ax.set_yticks(tick_positions)
        ax.set_yticklabels(tick_labels, fontsize=FONT_SIZE)

        ax.set_xlabel('Distance from anchor 1 (kb)', fontsize=FONT_SIZE, fontweight='bold')
        ax.set_ylabel('Distance from anchor 2 (kb)', fontsize=FONT_SIZE, fontweight='bold')

        # Get enrichment stat
        enrichment = result['statistics'].get('enrichment', 'N/A')
        if isinstance(enrichment, (int, float)):
            enrichment_str = f"{enrichment:.2f}x"
        else:
            enrichment_str = str(enrichment)

        ax.set_title(f'{category}\n(Enrichment: {enrichment_str})', fontsize=FONT_SIZE+2, fontweight='bold')

    # Remove empty subplots
    for idx in range(len(valid_results), len(axes)):
        fig.delaxes(axes[idx])

    # Add colorbar on the right side with sufficient space
    cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=FONT_SIZE-2)
    cbar.set_label('log10(contact frequency)', fontsize=FONT_SIZE, fontweight='bold', labelpad=20)

    plt.subplots_adjust(left=0.08, right=0.90, top=0.95, bottom=0.08, wspace=0.3, hspace=0.35)
    plt.savefig(os.path.join(output_dir, 'apa_comparison_heatmaps.png'), dpi=300, bbox_inches='tight')
    print(f"[visualizations][info] Saved APA comparison heatmap to {output_dir}/apa_comparison_heatmaps.png")
    plt.close()
