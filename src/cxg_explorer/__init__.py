# cellxgene-gene-explorer
# Copyright (c) 2026 TeraWeave Technologies. Licensed under the PolyForm Noncommercial License 1.0.0 (see LICENSE.md).
# Required Notice: Copyright (c) [YEAR] [YOUR NAME]
"""cellxgeneExplorer: single-gene analysis of CELLxGENE single-cell datasets, and cross-study summaries."""
__version__ = "1.0.0"


def explorer(out_dir="cellxgene_results", data_dir="cellxgene_data"):
    """Open the single-study panel in Jupyter."""
    from .notebook_agent import launch
    launch(out_dir, data_dir)


def summary_panel(results_folder="cellxgene_results"):
    """Open the cross-study summary panel in Jupyter."""
    from .notebook_summary import launch
    launch(results_folder)
