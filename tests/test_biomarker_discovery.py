import pytest
import numpy as np
import pandas as pd
from src.biomarker_discovery import DifferentialExpression, BiomarkerClassifier


def test_benjamini_hochberg():
    pvals = np.array([0.001, 0.01, 0.04, 0.2, 0.8])
    qvals = DifferentialExpression.benjamini_hochberg(pvals)
    assert len(qvals) == len(pvals)
    assert (qvals >= 0).all() and (qvals <= 1).all()
    # FDR düzeltmesi p-değerlerinden küçük olamaz
    assert (qvals >= pvals).all()


def test_differential_expression_analysis():
    np.random.seed(42)
    genes = ["GENE_A", "GENE_B", "GENE_C"]
    # Group A
    grp_a_data = np.random.normal(loc=5.0, scale=0.5, size=(20, 3))
    # Group B (GENE_A is lower, GENE_B is similar, GENE_C is higher)
    grp_b_data = np.random.normal(loc=2.0, scale=0.5, size=(20, 3))

    samples_a = [f"A_{i}" for i in range(20)]
    samples_b = [f"B_{i}" for i in range(20)]

    df = pd.DataFrame(
        np.vstack([grp_a_data, grp_b_data]),
        index=samples_a + samples_b,
        columns=genes,
    )

    res = DifferentialExpression.analyze(df, samples_a, samples_b, log2_fc_threshold=0.5, fdr_threshold=0.05)
    assert len(res) == 3
    assert "log2_fold_change" in res.columns
    assert "status" in res.columns


def test_biomarker_classifier():
    np.random.seed(42)
    n_samples = 40
    n_features = 10
    X = pd.DataFrame(
        np.random.randn(n_samples, n_features),
        columns=[f"Marker_{i}" for i in range(n_features)],
    )
    # y depends on Marker_0 and Marker_1
    y = pd.Series((X["Marker_0"] + X["Marker_1"] > 0).astype(int))

    clf = BiomarkerClassifier(n_estimators=30, random_state=42)
    importances = clf.fit_and_rank(X, y, top_k=3)
    assert len(importances) == n_features
    assert len(clf.selected_biomarkers) == 3

    cv_eval = clf.evaluate_cv(X, y, n_splits=3)
    assert "roc_auc" in cv_eval
    assert "accuracy" in cv_eval
    assert 0.0 <= cv_eval["roc_auc"] <= 1.0
