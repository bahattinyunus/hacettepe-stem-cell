"""
Multi-Omik Yolak Zenginleştirme (Pathway Enrichment & ORA) Analizörü.

Hipergeometrik dağılım ve Fisher's Exact Testi kullanarak diferansiyel
eksprese edilen genlerin kök hücre, pluripotensi, epigenetik ve onkolojik
yolaklardaki istatistiksel zenginleşmesini (Over-Representation Analysis) hesaplar.
"""

from typing import Dict, List, Optional, Tuple, Any
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


# Küratörlü Kök Hücre ve Rejeneratif Tıp Biyolojik Yolak Veritabanı
CURATED_PATHWAY_DATABASE = {
    "Pluripotensi ve Hücresel Yenilenme Ağı": [
        "POU5F1", "SOX2", "NANOG", "KLF4", "MYC", "ZFP42", "LIN28A", "SALL4", "ESRRB"
    ],
    "Wnt / Beta-Katenin Sinyalleşmesi & Kök Hücre Nişi": [
        "LGR5", "AXIN2", "CTNNB1", "MYC", "TCF7L2", "FZD7", "GSK3B", "DKK1"
    ],
    "Epitel-Mezenkim Geçişi (EMT) & Hücre Göçü": [
        "ACTA2", "COL1A1", "COL1A2", "SNAI1", "TWIST1", "CDH2", "FN1", "VIM", "FAP"
    ],
    "Glomerüler Filtrasyon & Podosit Bütünlüğü": [
        "NPHS1", "NPHS2", "WT1", "SYNPO", "ACTN4", "CD2AP", "PODXL"
    ],
    "Kemik İliği MSC İmmünomodülasyon & Sekretom": [
        "NT5E", "THY1", "ENG", "VEGFA", "HGF", "IGF1", "IL10", "TNFAIP6", "PTGS2"
    ],
    "DNA Hasar Yanıtı, Senescence & SASP": [
        "CDKN2A", "TP53", "ATM", "CHEK2", "IL6", "CXCL8", "SERPINE1", "GADD45A"
    ],
    "Kolorektal Kanser Kök Hücreliği & Kemorezistans": [
        "THY1", "PROM1", "CD44", "ALDH1A1", "LGR5", "ABCG2", "ABCB1", "EPCAM"
    ],
    "Hücresel Solunum & Mitokondriyal Metabolizma": [
        "GAPDH", "ACTB", "MT-CO1", "MT-ND4", "ATP5F1A", "COX4I1", "SOD2"
    ],
}


class PathwayEnrichment:
    """
    Gen setleri için Over-Representation Analysis (ORA) yürüten sınıf.
    """

    def __init__(self, pathway_db: Optional[Dict[str, List[str]]] = None):
        self.pathway_db = pathway_db or CURATED_PATHWAY_DATABASE

    @staticmethod
    def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
        """FDR düzeltmesi."""
        n = len(p_values)
        if n == 0:
            return np.array([])
        sorted_indices = np.argsort(p_values)
        sorted_p = p_values[sorted_indices]

        q_values = np.zeros(n)
        cummin = 1.0

        for i in range(n - 1, -1, -1):
            rank = i + 1
            q = sorted_p[i] * n / rank
            cummin = min(cummin, q)
            q_values[i] = min(cummin, 1.0)

        rev_indices = np.empty(n, dtype=int)
        rev_indices[sorted_indices] = np.arange(n)
        return q_values[rev_indices]

    def analyze(
        self,
        query_genes: List[str],
        background_universe_size: int = 500,
        min_overlap: int = 2,
    ) -> pd.DataFrame:
        """
        Sorgu gen listesinin veritabanındaki her bir yolakta zenginleşme derecesini hesaplar.

        Hipergeometrik test parametreleri:
        M = background_universe_size (toplam gen evreni)
        n = pathway_size (yolaktaki gen sayısı)
        N = len(query_genes) (sorgudaki gen sayısı)
        k = overlap_size (yolak ile sorgunun kesişimi)
        """
        results = []
        clean_query = set(g.strip().upper() for g in query_genes)
        N = len(clean_query)
        M = max(background_universe_size, 100)

        for pathway_name, genes in self.pathway_db.items():
            pathway_set = set(g.upper() for g in genes)
            overlap = clean_query.intersection(pathway_set)
            k = len(overlap)
            n = len(pathway_set)

            if k >= min_overlap:
                # P(X >= k) = sf(k - 1, M, n, N)
                p_val = hypergeom.sf(k - 1, M, n, N)
                rich_factor = k / n if n > 0 else 0.0

                # Odds ratio yaklaşık hesabı
                odds_ratio = (k * (M - n - N + k)) / ((n - k + 1e-4) * (N - k + 1e-4))

                results.append({
                    "pathway": pathway_name,
                    "pathway_size": n,
                    "overlap_count": k,
                    "overlap_genes": ", ".join(sorted(overlap)),
                    "rich_factor": rich_factor,
                    "odds_ratio": odds_ratio,
                    "p_value": float(p_val),
                })

        if not results:
            return pd.DataFrame(columns=[
                "pathway", "pathway_size", "overlap_count", "overlap_genes",
                "rich_factor", "odds_ratio", "p_value", "fdr_q_value", "-log10_p_value"
            ])

        df = pd.DataFrame(results)
        df["fdr_q_value"] = self.benjamini_hochberg(df["p_value"].values)
        df["-log10_p_value"] = -np.log10(np.maximum(df["p_value"].values, 1e-300))

        return df.sort_values(by="fdr_q_value").reset_index(drop=True)
