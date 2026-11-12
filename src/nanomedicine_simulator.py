"""
Biyomimetik Nanotıp ve Kanser Kök Hücresi Hedefleme Simülatörü.

Araş. Gör. Dr. Özlem ALTUNDAĞ ERDOĞAN'ın çalıştay sunumunda detaylandırılan
hücre zarıyla kaplı (cell membrane-coated) akıllı nanotaşıyıcıların tümör
mikroçevresine nüfuz etme, aktif THY1/CD90 bağlanma ve hücre içi ilaç salım
farmakokinetik/farmakodinamik (PK/PD) 2-kompartmanlı simülasyonu.
"""

from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd


class BiomimeticNanomedicine:
    """
    Biyomimetik nanotaşıyıcı ve serbest kemoterapi ajanının
    tümör ve kanser kök hücresi hedefleme kinetiğini simüle eden sınıf.
    """

    def __init__(
        self,
        drug_name: str = "5-Fluorouracil (5-FU)",
        membrane_coating: str = "Mezenkimal Kök Hücre (MSC) Zarı",
        targeting_ligand: str = "Anti-THY1 (CD90) Aptameri",
    ):
        self.drug_name = drug_name
        self.membrane_coating = membrane_coating
        self.targeting_ligand = targeting_ligand

    def simulate_pharmacokinetics(
        self,
        dose_mg_kg: float = 20.0,
        time_hours: float = 48.0,
        n_points: int = 100,
        stealth_coating_efficiency: float = 0.90,
    ) -> Dict[str, Any]:
        """
        2-kompartmanlı farmakokinetik model:
        C_plasma(t) = A * exp(-alpha * t) + B * exp(-beta * t)
        Tümör dokusu konsantrasyonu: dC_tumor/dt = k_in * C_plasma - k_out * C_tumor
        """
        t = np.linspace(0, time_hours, n_points)

        # Serbest ilaç: Hızlı klerens (yarı ömür ~ 0.5 - 2 saat)
        alpha_free = 1.4
        beta_free = 0.35
        plasma_free = dose_mg_kg * (0.8 * np.exp(-alpha_free * t) + 0.2 * np.exp(-beta_free * t))

        # Biyomimetik taşıyıcı: RES kaçışı sayesinde uzun dolaşım (yarı ömür ~ 12 - 24 saat)
        clearance_reduction = np.clip(stealth_coating_efficiency, 0.2, 0.98)
        alpha_nano = alpha_free * (1.0 - clearance_reduction * 0.75)
        beta_nano = beta_free * (1.0 - clearance_reduction * 0.85)
        plasma_nano = dose_mg_kg * (0.6 * np.exp(-alpha_nano * t) + 0.4 * np.exp(-beta_nano * t))

        # Tümör akümülasyonu (EPR etkisi ve aktif kargo tutunması)
        k_in_free = 0.08
        k_out_free = 0.22
        k_in_nano = 0.45 * (1.0 + stealth_coating_efficiency * 0.5)
        k_out_nano = 0.05  # Nanopartiküller tümörde hapsolur

        tumor_free = np.zeros_like(t)
        tumor_nano = np.zeros_like(t)

        dt = t[1] - t[0]
        for i in range(1, len(t)):
            tumor_free[i] = tumor_free[i-1] + (k_in_free * plasma_free[i-1] - k_out_free * tumor_free[i-1]) * dt
            tumor_nano[i] = tumor_nano[i-1] + (k_in_nano * plasma_nano[i-1] - k_out_nano * tumor_nano[i-1]) * dt

        # CSC içi pH duyarlı salım (asidik lizozomda serbestleşen kargo)
        lysosomal_release_efficiency = 0.88
        csc_intracellular_drug = tumor_nano * lysosomal_release_efficiency

        # Trapezoidal numerical integration (NumPy 2.0+ uyumlu)
        trapz_func = getattr(np, "trapezoid", getattr(np, "trapz", None))
        auc_free_tumor = float(trapz_func(tumor_free, t))
        auc_nano_tumor = float(trapz_func(tumor_nano, t))
        tumor_targeting_fold_increase = auc_nano_tumor / (auc_free_tumor + 1e-6)

        return {
            "time_hours": t,
            "plasma_free": plasma_free,
            "plasma_nano": plasma_nano,
            "tumor_free": tumor_free,
            "tumor_nano": tumor_nano,
            "csc_intracellular_drug": csc_intracellular_drug,
            "auc_free_tumor": auc_free_tumor,
            "auc_nano_tumor": auc_nano_tumor,
            "tumor_targeting_fold_increase": tumor_targeting_fold_increase,
        }

    def simulate_dose_response(
        self,
        concentrations: Optional[np.ndarray] = None,
    ) -> pd.DataFrame:
        """
        Kanser kök hücreleri (CSC) ve normal somatik hücreler üzerinde
        serbest kemoterapi vs biyomimetik formülasyonun hücre sağkalım eğrisi (Hill denklemi).
        """
        if concentrations is None:
            concentrations = np.logspace(-2, 3, 50)  # 0.01 - 1000 uM

        # Serbest kemoterapide CSC'ler kemorezistandır (yüksek IC50)
        ic50_csc_free = 145.0  # uM
        # Biyomimetik hedefli partikülde CSC IC50 dramatik düşer
        ic50_csc_nano = 6.5    # uM
        # Normal somatik hücreler
        ic50_normal_free = 18.0  # Toksik
        ic50_normal_nano = 280.0 # Nanopartikül kalkanı sayesinde somatik hücreye toksisite düşük

        hill_coef = 1.6

        viability_csc_free = 100.0 / (1.0 + (concentrations / ic50_csc_free) ** hill_coef)
        viability_csc_nano = 100.0 / (1.0 + (concentrations / ic50_csc_nano) ** hill_coef)
        viability_normal_free = 100.0 / (1.0 + (concentrations / ic50_normal_free) ** hill_coef)
        viability_normal_nano = 100.0 / (1.0 + (concentrations / ic50_normal_nano) ** hill_coef)

        therapeutic_index_free = ic50_normal_free / ic50_csc_free
        therapeutic_index_nano = ic50_normal_nano / ic50_csc_nano

        return pd.DataFrame({
            "concentration_uM": concentrations,
            "CSC_Free_Drug": viability_csc_free,
            "CSC_Biomimetic_Nano": viability_csc_nano,
            "Normal_Free_Drug": viability_normal_free,
            "Normal_Biomimetic_Nano": viability_normal_nano,
            "TI_Free": therapeutic_index_free,
            "TI_Nano": therapeutic_index_nano,
        })
