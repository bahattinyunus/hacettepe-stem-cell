import pytest
from src.secretome_network import SecretomeNetwork, PARACRINE_DATABASE


def test_secretome_network_dataframe():
    sec_net = SecretomeNetwork()
    df = sec_net.get_dataframe()
    assert len(df) == len(PARACRINE_DATABASE)
    assert "ligand" in df.columns
    assert "receptor" in df.columns
    assert "target_cell" in df.columns


def test_secretome_filter_by_organ():
    sec_net = SecretomeNetwork()
    kidney_df = sec_net.filter_by_organ("Böbrek")
    assert len(kidney_df) > 0
    assert all("Böbrek" in org for org in kidney_df["organ"])


def test_simulate_podocyte_rescue():
    sec_net = SecretomeNetwork()
    sim = sec_net.simulate_podocyte_rescue(msc_dose=1.5, exosome_enriched=True, initial_damage_score=0.8)
    assert sim["final_damage"] < sim["initial_damage"]
    assert sim["nephrin_recovery"] > 0.5
    assert sim["proteinuria_reduction_pct"] > 0
    assert len(sim["timeline"]["days"]) == 29
