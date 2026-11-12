import subprocess
import sys
import pytest


def test_cli_help():
    cmd = [sys.executable, "-m", "src.cli", "--help"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    assert "simulate-reprogramming" in res.stdout
    assert "calc-stemness" in res.stdout
    assert "epigenetic-clock" in res.stdout


def test_cli_epigenetic_clock():
    cmd = [sys.executable, "-m", "src.cli", "epigenetic-clock", "--age", "60", "--cycles", "4"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    assert "Kronolojik Yas: 60.0" in res.stdout or "Kronolojik Yaş: 60.0" in res.stdout
    assert "Rejuvenasyon Sonras" in res.stdout


def test_cli_pathway_enrichment():
    cmd = [sys.executable, "-m", "src.cli", "pathway-enrichment", "--genes", "POU5F1,SOX2,NANOG"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    assert res.returncode == 0
    assert "Pluripotensi" in res.stdout
