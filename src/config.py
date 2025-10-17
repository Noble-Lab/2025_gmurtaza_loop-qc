"""
Configuration file handling for loop analysis pipeline.

Loads and validates YAML configuration files.
"""

import yaml
import os
import sys


def load_config(config_path):
    """
    Load and validate configuration from YAML file.

    Parameters:
    -----------
    config_path : str
        Path to YAML configuration file

    Returns:
    --------
    dict
        Configuration dictionary
    """
    if not os.path.exists(config_path):
        print(f"[config][error] Configuration file not found: {config_path}")
        sys.exit(1)

    try:
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
    except yaml.YAMLError as e:
        print(f"[config][error] Failed to parse YAML configuration: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"[config][error] Error reading configuration file: {e}")
        sys.exit(1)

    # Validate required sections
    required_sections = ['input', 'output', 'loops']
    for section in required_sections:
        if section not in config:
            print(f"[config][error] Missing required section: {section}")
            sys.exit(1)

    # Validate input section
    if 'bedpe' not in config['input']:
        print(f"[config][error] Missing 'bedpe' in input section")
        sys.exit(1)

    if 'directory' not in config['output']:
        print(f"[config][error] Missing 'directory' in output section")
        sys.exit(1)

    # Check BEDPE file exists
    bedpe_path = config['input']['bedpe']
    if not os.path.exists(bedpe_path):
        print(f"[config][error] BEDPE file not found: {bedpe_path}")
        sys.exit(1)

    print(f"[config][info] Configuration loaded successfully from: {config_path}")

    return config


def get_motif_config(config):
    """
    Extract and validate motif analysis configuration.

    Parameters:
    -----------
    config : dict
        Configuration dictionary

    Returns:
    --------
    dict or None
        Motif analysis config with 'fasta' and 'motif' keys, or None if not configured
    """
    motif_config = config.get('motif_analysis')

    if motif_config is None:
        return None

    # Check if both fasta and motif are specified
    fasta = motif_config.get('fasta')
    motif = motif_config.get('motif')

    if fasta is None or motif is None:
        return None

    # Validate files exist
    if not os.path.exists(fasta):
        print(f"[config][error] Reference FASTA file not found: {fasta}")
        sys.exit(1)

    if not os.path.exists(motif):
        print(f"[config][error] Motif file not found: {motif}")
        sys.exit(1)

    print(f"[config][info] Motif analysis configured")

    return {
        'fasta': fasta,
        'motif': motif,
        'fdr_threshold': motif_config.get('fdr_threshold', 0.05),
        'fimo_threshold': motif_config.get('fimo_threshold', '1e-4'),
        'optimize_convergent': motif_config.get('optimize_convergent', False)
    }


def get_apa_config(config):
    """
    Extract and validate APA analysis configuration.

    Parameters:
    -----------
    config : dict
        Configuration dictionary

    Returns:
    --------
    dict or None
        APA analysis config with 'hic_file', 'window_size', and 'resolution' keys, or None if not configured
    """
    apa_config = config.get('apa_analysis')

    if apa_config is None:
        return None

    # Check if hic_file is specified
    hic_file = apa_config.get('hic_file')

    if hic_file is None:
        return None

    # Validate Hi-C file exists
    if not os.path.exists(hic_file):
        print(f"[config][error] Hi-C file not found: {hic_file}")
        sys.exit(1)

    print(f"[config][info] APA analysis configured")

    return {
        'hic_file': hic_file,
        'window_size': apa_config.get('window_size', 100000),
        'resolution': apa_config.get('resolution', 10000),
        'nproc': apa_config.get('nproc', 1),
        'min_loop_length': apa_config.get('min_loop_length', None)
    }
