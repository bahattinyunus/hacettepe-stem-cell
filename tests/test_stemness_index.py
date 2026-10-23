import pytest
import numpy as np
import pandas as pd
from src.stemness_index import calculate_stemness_index, StemnessModel
from src.scrna_pipeline import generate_synthetic_stem_cell_data, ScRNAPipeline


def test_stemness_index_calculation():
    counts_df, metadata_df = generate_synthetic_stem_cell_data(n_cells=80, seed=42)
    pipeline = ScRNAPipeline()
    log1p_df = pipeline.preprocess(counts_df)

    stemness = calculate_stemness_index(log1p_df)
    assert len(stemness) == len(log1p_df)
    assert stemness.min() >= 0.0
    assert stemness.max() <= 1.0

    # Pluripotent ESC'lerin kök hücrelilik indeksi genellikle somatik fibroblastlardan yüksektir
    esc_cells = metadata_df[metadata_df["cell_type"] == "Pluripotent ESC"].index.intersection(stemness.index)
    fib_cells = metadata_df[metadata_df["cell_type"] == "Somatik Fibroblast"].index.intersection(stemness.index)

    if len(esc_cells) > 0 and len(fib_cells) > 0:
        mean_esc = stemness.loc[esc_cells].mean()
        mean_fib = stemness.loc[fib_cells].mean()
        assert mean_esc > mean_fib
