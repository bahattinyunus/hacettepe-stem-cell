import numpy as np
import pytest
from src.waddington import WaddingtonLandscape, simulate_reprogramming


def test_waddington_potential_shape_and_values():
    landscape = WaddingtonLandscape()
    x = np.array([-1.5, 0.0, 1.5])
    y = np.array([0.0, 0.0, 0.0])
    pot = landscape.potential(x, y, yamanaka_factor_strength=0.0)
    assert len(pot) == 3
    assert not np.isnan(pot).any()


def test_waddington_gradient():
    landscape = WaddingtonLandscape()
    gx, gy = landscape.gradient(0.5, 0.2, yamanaka_factor_strength=0.5)
    assert isinstance(gx, float)
    assert isinstance(gy, float)
    assert not np.isnan(gx) and not np.isnan(gy)


def test_waddington_surface_grid():
    landscape = WaddingtonLandscape()
    xx, yy, zz = landscape.generate_surface_grid(grid_points=20)
    assert xx.shape == (20, 20)
    assert yy.shape == (20, 20)
    assert zz.shape == (20, 20)


def test_simulate_reprogramming():
    result = simulate_reprogramming(start_pos=(1.3, 0.6), yamanaka_strength=0.9, steps=50, seed=42)
    assert "trajectory" in result
    assert "energy" in result
    assert "states" in result
    assert result["trajectory"].shape == (50, 2)
    assert len(result["states"]) == 50
    assert isinstance(result["success"], bool)
