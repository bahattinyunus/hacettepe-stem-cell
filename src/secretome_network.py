"""
Mezenkimal Kök Hücre (MSC) Sekretom ve Parakrin Mikroçevre Ağı Modülü.

Kemik iliği MSC'lerinin hasarlı böbrek podositleri ve fibrotik mikroçevrede
salgıladığı parakrin faktörlerin (büyüme faktörleri, antienflamatuvar sitokinler, eksozomlar)
hedef hücre reseptörleriyle etkileşimini ve rejeneratif yanıtını modeller.
"""

from typing import Dict, List, Optional, Tuple, Any
import pandas as pd
import numpy as np


# Hacettepe Kök Hücre Çalıştayı II. Oturumu temel alınarak derlenen parakrin veritabanı
PARACRINE_DATABASE = [
    {
        "ligand": "VEGFA",
        "ligand_source": "Kemik Iligi MSC",
        "receptor": "KDR (VEGFR2)",
        "target_cell": "Hasarli Podosit / Endotel",
        "pathway": "Anjiyogenez & Endotelyal Koruma",
        "action": "Podosit ayaksı çıkıntılarının (foot processes) stabilizasyonu",
        "weight": 0.95,
        "organ": "Böbrek Glomerülü",
    },
    {
        "ligand": "HGF",
        "ligand_source": "Kemik Iligi MSC",
        "receptor": "MET (c-Met)",
        "target_cell": "Hasarli Podosit",
        "pathway": "Anti-Apoptoz & Doku Rejenerasyonu",
        "action": "Kaspaz-3 inhibisyonu ve hücre sağkalımı",
        "weight": 0.90,
        "organ": "Böbrek Glomerülü",
    },
    {
        "ligand": "IGF1",
        "ligand_source": "Kemik Iligi MSC",
        "receptor": "IGF1R",
        "target_cell": "Podosit",
        "pathway": "PI3K/Akt Pro-Survival",
        "action": "NPHS1 (Nefrin) ekspresyonunun restorasyonu",
        "weight": 0.85,
        "organ": "Böbrek Glomerülü",
    },
    {
        "ligand": "IL10",
        "ligand_source": "Kemik Iligi MSC / Eksozom",
        "receptor": "IL10RA",
        "target_cell": "Makrofaj / T Hücresi",
        "pathway": "İmmünomodülasyon",
        "action": "M1 pro-inflamatuvar makrofajların M2 iyileştirici fenotipe dönüşümü",
        "weight": 0.88,
        "organ": "İnflamatuvar Mikroçevre",
    },
    {
        "ligand": "TNFAIP6 (TSG-6)",
        "ligand_source": "Kemik Iligi MSC",
        "receptor": "CD44",
        "target_cell": "İnflamatuvar İnfiltratlar",
        "pathway": "Anti-Enflamasyon",
        "action": "NF-kB yolağının baskılanması ve sitokin fırtınasının durdurulması",
        "weight": 0.92,
        "organ": "Kemik İliği / Akciğer / Böbrek",
    },
    {
        "ligand": "FGF2 (bFGF)",
        "ligand_source": "Kemik Iligi MSC",
        "receptor": "FGFR1",
        "target_cell": "Fibroblast / Alveoler Hücre",
        "pathway": "Fibrozis Modülasyonu",
        "action": "Pulmoner fibroziste aşırı kollajen birikiminin sınırlandırılması",
        "weight": 0.82,
        "organ": "Akciğer Stroma",
    },
    {
        "ligand": "ANGPT1 (Angiopoietin-1)",
        "ligand_source": "Kemik Iligi MSC",
        "receptor": "TEK (Tie2)",
        "target_cell": "Glomerüler Endotel",
        "pathway": "Vasküler Permeabilite Bariyeri",
        "action": "Proteinürinin önlenmesi ve filtrasyon bariyer bütünlüğü",
        "weight": 0.86,
        "organ": "Böbrek Glomerülü",
    },
]


class SecretomeNetwork:
    """
    MSC parakrin sinyalleşmesini ve doku kurtarma (rescue) skorunu hesaplayan sınıf.
    """

    def __init__(self, db: Optional[List[Dict[str, Any]]] = None):
        self.db = db or PARACRINE_DATABASE
        self.df = pd.DataFrame(self.db)

    def get_dataframe(self) -> pd.DataFrame:
        """Etkileşim tablosunu DataFrame olarak döndürür."""
        return self.df.copy()

    def filter_by_organ(self, organ_keyword: str) -> pd.DataFrame:
        """Organ veya doku tipine göre sinyalleşme akslarını filtreler."""
        mask = self.df["organ"].str.contains(organ_keyword, case=False, na=False)
        return self.df[mask]

    def simulate_podocyte_rescue(
        self,
        msc_dose: float = 1.0,
        exosome_enriched: bool = True,
        initial_damage_score: float = 0.85,
    ) -> Dict[str, Any]:
        """
        MSC uygulaması ve eksozom tedavisi sonrasında podosit mikroçevresindeki
        iyileşme oranını, proteinüri düşüşünü ve nefrin/podosin geri kazanımını hesaplar.

        Args:
            msc_dose: MSC hücre dozu faktörü (0.1 - 2.5)
            exosome_enriched: Eksozom saflaştırılmış formülasyon varlığı
            initial_damage_score: Başlangıçtaki glomerüler hasar derecesi [0 - 1]

        Returns:
            Simülasyon sonuçları sözlüğü.
        """
        exosome_multiplier = 1.35 if exosome_enriched else 1.0
        rescue_power = np.clip(0.35 * msc_dose * exosome_multiplier, 0.05, 0.95)

        residual_damage = max(0.05, initial_damage_score * (1.0 - rescue_power))
        nephrin_recovery = min(1.0, (1.0 - residual_damage) * 1.15)
        podocin_recovery = min(1.0, (1.0 - residual_damage) * 1.10)
        proteinuria_reduction_pct = (initial_damage_score - residual_damage) / initial_damage_score * 100

        # Zaman serisi simülasyonu (28 günlük izlem)
        days = np.arange(0, 29)
        # Lojistik iyileşme eğrisi
        k = 0.22 * msc_dose * exosome_multiplier
        t0 = 7.0
        time_curve_damage = residual_damage + (initial_damage_score - residual_damage) / (1 + np.exp(k * (days - t0)))
        time_curve_nephrin = 0.2 + (nephrin_recovery - 0.2) / (1 + np.exp(-k * (days - t0)))

        return {
            "msc_dose": msc_dose,
            "exosome_enriched": exosome_enriched,
            "initial_damage": initial_damage_score,
            "final_damage": residual_damage,
            "nephrin_recovery": nephrin_recovery,
            "podocin_recovery": podocin_recovery,
            "proteinuria_reduction_pct": proteinuria_reduction_pct,
            "timeline": {
                "days": days,
                "damage_curve": time_curve_damage,
                "nephrin_curve": time_curve_nephrin,
            },
        }
