import pytest
import pandas as pd
import numpy as np
from src.scrna_pipeline import (
    generate_synthetic_stem_cell_data,
    ScRNAPipeline,
)


def test_generate_synthetic_stem_cell_data():
    counts_df, metadata_df = generate_synthetic_stem_cell_data(n_cells=100, n_genes=35, seed=42)
    assert counts_df.shape[0] == 100
    assert counts_df.shape[1] == 35
    assert metadata_df.shape[0] == 100
    assert "cell_type" in metadata_df.columns
    assert "total_counts" in metadata_df.columns


def test_scrna_pipeline_workflow():
    counts_df, metadata_df = generate_synthetic_stem_cell_data(n_cells=80, n_genes=30, seed=42)
    pipeline = ScRNAPipeline()

    # Preprocess
    log1p_df = pipeline.preprocess(counts_df, min_counts=10, min_genes=5)
    assert log1p_df.shape[0] <= 80
    assert (log1p_df >= 0).all().all()

    # HVG
    hvgs = pipeline.find_highly_variable_genes(top_n=15)
    assert len(hvgs) == 15

    # PCA
    pca_df = pipeline.run_pca(n_components=3)
    assert pca_df.shape[1] == 3
    assert pca_df.shape[0] == log1p_df.shape[0]

    # Clusters
    clusters = pipeline.cluster_cells(n_clusters=3)
    assert len(clusters) == len(pca_df)
    assert set(clusters.unique()).issubset({0, 1, 2})

    # Pseudotime
    pseudotime = pipeline.compute_pseudotime()
    assert len(pseudotime) == len(pca_df)
    assert pseudotime.min() >= 0.0
    assert pseudotime.max() <= 1.0


def test_differential_expression_gene():
    counts_df, metadata_df = generate_synthetic_stem_cell_data(n_cells=120, seed=42)
    pipeline = ScRNAPipeline()
    pipeline.preprocess(counts_df)

    esc_cells = metadata_df[metadata_df["cell_type"] == "Pluripotent ESC"].index.tolist()
    fibro_cells = metadata_df[metadata_df["cell_type"] == "Somatik Fibroblast"].index.tolist()

    res = pipeline.differential_expression_gene("POU5F1", esc_cells, fibro_cells)
    assert res["gene"] == "POU5F1"
    assert res["log2_fold_change"] > 0
    assert "p_value" in res
