import pytest
from src.pathway_enrichment import PathwayEnrichment


def test_pathway_enrichment_pluripotency():
    ora = PathwayEnrichment()
    query = ["POU5F1", "SOX2", "NANOG", "KLF4"]
    df = ora.analyze(query, background_universe_size=300, min_overlap=2)

    assert not df.empty
    assert "Pluripotensi ve Hücresel Yenilenme Ağı" in df["pathway"].values
    top_row = df[df["pathway"] == "Pluripotensi ve Hücresel Yenilenme Ağı"].iloc[0]
    assert top_row["overlap_count"] == 4
    assert top_row["p_value"] < 0.05
    assert top_row["fdr_q_value"] < 0.05


def test_pathway_enrichment_empty():
    ora = PathwayEnrichment()
    df = ora.analyze(["NON_EXISTENT_GENE_123"], min_overlap=2)
    assert df.empty
