"""
Unit tests for the In Silico Whole-Cell Digital Twin & Perturbation Engine.
"""

import pytest
from src.digital_twin import WholeCellDigitalTwin, PerturbationResult


def test_digital_twin_initialization():
    twin = WholeCellDigitalTwin()
    assert len(twin.genes) >= 15
    assert "POU5F1" in twin.genes
    assert "SETD2" in twin.genes
    assert "NPHS1" in twin.genes


def test_setd2_epigenetic_barrier_knockout():
    twin = WholeCellDigitalTwin()
    res = twin.simulate_perturbation("SETD2", perturbation_type="CRISPR_KO", intensity=0.95)
    
    assert isinstance(res, PerturbationResult)
    assert res.target == "SETD2"
    assert res.reprogramming_efficiency_multiplier > 10.0
    assert res.fold_changes["SETD2"] < 0
    assert res.fold_changes["POU5F1"] > 0
    assert "H3K36me3" in res.biomarker_alert


def test_gyy4137_pulmonary_fibrosis_rescue():
    twin = WholeCellDigitalTwin()
    res = twin.simulate_perturbation("GYY4137", perturbation_type="TROPHIC_RESCUE", intensity=0.85)
    
    assert "Akciğer" in res.phenotypic_transition
    assert res.fold_changes["IL1B"] < 0
    assert res.viability_percentage > 70.0


def test_niclosamide_csc_inhibition():
    twin = WholeCellDigitalTwin()
    res = twin.simulate_perturbation("Niclosamide", perturbation_type="INHIBITION", intensity=0.90)
    
    assert "Kanser Kök Hücresi" in res.phenotypic_transition
    assert res.fold_changes["ABCG1"] < 0
    assert res.fold_changes["MMP2"] < 0


def test_dataframe_conversion():
    twin = WholeCellDigitalTwin()
    res = twin.simulate_perturbation("PODOCYTE_RESCUE", perturbation_type="TROPHIC_RESCUE")
    df = twin.to_dataframe(res)
    
    assert len(df) == len(twin.genes)
    assert "Log2FC" in df.columns
    assert "Significant" in df.columns
