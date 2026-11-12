"""
Epigenetik Saat (Horvath DNA Metilasyonu) ve In Vivo Rejuvenasyon Modülü.

Steve Horvath'ın (Genome Biology, 2013) pan-dokusal epigenetik saati ve
yeni nesil hücresel gençleşme (pulsed OSKM rejuvenation) dinamiklerinin
matematiksel ve biyoinformatik modellemesi.
"""

from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd


# Horvath Saati için kritik CpG ve epigenetik belirteç gen lokusları
KEY_EPIGENETIC_CPG_LOCI = {
    "ELOVL2": {"coef": 0.085, "type": "hypermethylated_with_age", "baseline_young": 0.15, "baseline_old": 0.82},
    "FHL2": {"coef": 0.078, "type": "hypermethylated_with_age", "baseline_young": 0.10, "baseline_old": 0.75},
    "PENK": {"coef": 0.062, "type": "hypermethylated_with_age", "baseline_young": 0.12, "baseline_old": 0.68},
    "GRIA2": {"coef": 0.055, "type": "hypermethylated_with_age", "baseline_young": 0.20, "baseline_old": 0.72},
    "KLF14": {"coef": 0.048, "type": "hypermethylated_with_age", "baseline_young": 0.08, "baseline_old": 0.60},
    "CDKN2A_p16": {"coef": 0.092, "type": "hypermethylated_with_age", "baseline_young": 0.05, "baseline_old": 0.85},
    # Yaşla birlikte hipometillenen (koruyucu etkisi azalan) lokuslar:
    "TERT_promoter": {"coef": -0.065, "type": "hypomethylated_with_age", "baseline_young": 0.78, "baseline_old": 0.25},
    "DNMT3A": {"coef": -0.058, "type": "hypomethylated_with_age", "baseline_young": 0.82, "baseline_old": 0.35},
    "TET2": {"coef": -0.052, "type": "hypomethylated_with_age", "baseline_young": 0.75, "baseline_old": 0.30},
}


class HorvathEpigeneticClock:
    """
    DNA metilasyon beta değerlerinden biyolojik yaşı (DNAmAge) ve
    yaşlanma ivmesini (Age Acceleration) hesaplayan model.
    """

    def __init__(self, loci_dict: Optional[Dict[str, Dict[str, Any]]] = None):
        self.loci = loci_dict or KEY_EPIGENETIC_CPG_LOCI
        self.intercept = 20.0

    def generate_synthetic_methylation_profile(
        self,
        chronological_age: float,
        lifestyle_stress_factor: float = 0.0,
        seed: Optional[int] = 42,
    ) -> Dict[str, float]:
        """
        Belirli bir kronolojik yaş ve stres faktörüne göre sentetik CpG beta değerleri üretir.
        """
        if seed is not None:
            np.random.seed(seed)

        # 0 - 100 yaş arası yaş fraksiyonu
        age_frac = np.clip(chronological_age / 100.0, 0.0, 1.0)
        profile = {}

        for gene, info in self.loci.items():
            noise = np.random.normal(0, 0.03)
            stress_impact = lifestyle_stress_factor * 0.05

            if info["type"] == "hypermethylated_with_age":
                # Yaşlandıkça metilasyon artar
                beta = info["baseline_young"] + age_frac * (info["baseline_old"] - info["baseline_young"]) + stress_impact + noise
            else:
                # Yaşlandıkça metilasyon düşer
                beta = info["baseline_young"] - age_frac * (info["baseline_young"] - info["baseline_old"]) - stress_impact + noise

            profile[gene] = float(np.clip(beta, 0.01, 0.99))

        return profile

    def predict_biological_age(self, beta_values: Dict[str, float]) -> float:
        """
        Horvath algoritması benzeri ağırlıklı toplam ile DNA metilasyon yaşını (DNAmAge) hesaplar.
        """
        weighted_sum = 0.0
        total_weight = 0.0

        for gene, info in self.loci.items():
            if gene in beta_values:
                beta = beta_values[gene]
                coef = info["coef"]
                if coef > 0:
                    # Hipermetilasyon: Yüksek beta -> Yaşlı
                    contrib = (beta - info["baseline_young"]) / (info["baseline_old"] - info["baseline_young"] + 1e-6)
                else:
                    # Hipometilasyon: Düşük beta -> Yaşlı
                    contrib = (info["baseline_young"] - beta) / (info["baseline_young"] - info["baseline_old"] + 1e-6)

                weighted_sum += abs(coef) * contrib * 100.0
                total_weight += abs(coef)

        if total_weight == 0:
            return self.intercept

        dnam_age = weighted_sum / total_weight
        return float(np.clip(dnam_age, 0.0, 120.0))

    def simulate_pulsed_rejuvenation(
        self,
        chronological_age: float = 65.0,
        cycles: int = 6,
        pulse_strength: float = 0.85,
        preservation_safety_threshold: float = 0.90,
    ) -> Dict[str, Any]:
        """
        In vivo geçici (pulsatile) Yamanaka faktörü ekspresyonunun
        epigenetik saati sıfırlama ve hücresel kimliği koruma simülasyonu.

        Args:
            chronological_age: Bireyin mevcut kronolojik yaşı
            cycles: 2 gün açık / 5 gün kapalı protokol döngü sayısı
            pulse_strength: Yamanaka OSKM dozu (0.1 - 1.0)
            preservation_safety_threshold: Somatik kimliğin korunma güvenliği eşiği
        """
        initial_profile = self.generate_synthetic_methylation_profile(chronological_age)
        initial_bio_age = self.predict_biological_age(initial_profile)

        timeline_cycles = np.arange(0, cycles + 1)
        bio_age_curve = [initial_bio_age]
        identity_loss_risk = [0.0]

        current_profile = initial_profile.copy()

        for c in range(1, cycles + 1):
            # Her döngüde epigenetik gürültü temizlenir ve beta değerleri gençleşme yönünde evrilir
            reversal_rate = 0.12 * pulse_strength
            for gene, info in self.loci.items():
                target_young = info["baseline_young"]
                curr_val = current_profile[gene]
                # Gençlik zeminine doğru gerileme
                current_profile[gene] = curr_val + reversal_rate * (target_young - curr_val)

            new_bio_age = self.predict_biological_age(current_profile)
            bio_age_curve.append(new_bio_age)

            # Teratom / pluripotensiye aşırı kaçış riski (pulse strength çok yüksekse ve döngü uzarsa artar)
            risk = np.clip((c * (pulse_strength ** 2) * 0.035) if pulse_strength > 0.80 else 0.01 * c, 0.0, 1.0)
            identity_loss_risk.append(risk)

        final_bio_age = bio_age_curve[-1]
        years_rejuvenated = max(0.0, initial_bio_age - final_bio_age)
        safety_status = "Güvenli (Somatik Kimlik Korundu)" if identity_loss_risk[-1] < (1.0 - preservation_safety_threshold) else "Riskli (Dediferansiasyon / Teratom Tehlikesi)"

        return {
            "chronological_age": chronological_age,
            "initial_biological_age": initial_bio_age,
            "final_biological_age": final_bio_age,
            "years_rejuvenated": years_rejuvenated,
            "cycles": cycles,
            "timeline": {
                "cycles": timeline_cycles,
                "bio_age_curve": np.array(bio_age_curve),
                "identity_loss_risk": np.array(identity_loss_risk),
            },
            "safety_status": safety_status,
            "rejuvenated_profile": current_profile,
        }
