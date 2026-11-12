import pytest
from src.nanomedicine_simulator import BiomimeticNanomedicine


def test_nanomedicine_pharmacokinetics():
    nano = BiomimeticNanomedicine()
    pk = nano.simulate_pharmacokinetics(dose_mg_kg=25.0, time_hours=36.0, n_points=50)

    assert "time_hours" in pk
    assert "tumor_nano" in pk
    assert "tumor_free" in pk
    assert pk["auc_nano_tumor"] > pk["auc_free_tumor"]
    assert pk["tumor_targeting_fold_increase"] > 1.5


def test_nanomedicine_dose_response():
    nano = BiomimeticNanomedicine()
    dr = nano.simulate_dose_response()

    assert not dr.empty
    assert "CSC_Free_Drug" in dr.columns
    assert "CSC_Biomimetic_Nano" in dr.columns
    assert dr["TI_Nano"].iloc[0] > dr["TI_Free"].iloc[0]
