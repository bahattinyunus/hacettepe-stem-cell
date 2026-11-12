import pytest
from src.epigenetic_clock import HorvathEpigeneticClock


def test_epigenetic_clock_prediction():
    clock = HorvathEpigeneticClock()
    young_profile = clock.generate_synthetic_methylation_profile(chronological_age=25.0)
    old_profile = clock.generate_synthetic_methylation_profile(chronological_age=75.0)

    age_young = clock.predict_biological_age(young_profile)
    age_old = clock.predict_biological_age(old_profile)

    assert age_young < age_old
    assert 15.0 <= age_young <= 40.0
    assert 55.0 <= age_old <= 95.0


def test_simulate_pulsed_rejuvenation():
    clock = HorvathEpigeneticClock()
    res = clock.simulate_pulsed_rejuvenation(chronological_age=68.0, cycles=5, pulse_strength=0.75)

    assert res["chronological_age"] == 68.0
    assert res["final_biological_age"] < res["initial_biological_age"]
    assert res["years_rejuvenated"] > 0
    assert len(res["timeline"]["cycles"]) == 6
    assert "safety_status" in res
