"""
Kök Hücrelilik İndeksi (Stemness Index - mRNAsI) Hesaplama Modülü.

Malta et al. (Cell, 2018) tarafından geliştirilen ve kök hücre/pluripotensi
imzasını (stemness signature) kullanarak hücrelerin farklılaşmamışlık
derecesini [0, 1] skalasında puanlayan algoritmanın Python uygulaması.
"""

from typing import Dict, List, Optional, Union
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


# Pluripotensi ve kök hücrelilik ağırlık vektörü (referans genler)
DEFAULT_STEMNESS_WEIGHTS: Dict[str, float] = {
    "POU5F1": 0.85,
    "SOX2": 0.80,
    "NANOG": 0.82,
    "KLF4": 0.65,
    "MYC": 0.55,
    "ZFP42": 0.70,
    "LIN28A": 0.75,
    "THY1": 0.50,
    "PROM1": 0.60,
    "CD44": 0.45,
    "ALDH1A1": 0.48,
    "LGR5": 0.62,
    "EPCAM": 0.40,
    # Farklılaşma yönündeki genler negatif katsayı alır
    "COL1A1": -0.70,
    "COL1A2": -0.65,
    "ACTA2": -0.60,
    "NPHS1": -0.50,
    "PTPRC": -0.55,
}


class StemnessModel:
    """
    Ekspresyon profillerinden kök hücrelilik indeksini (mRNAsI) hesaplayan sınıf.
    """

    def __init__(self, weights: Optional[Dict[str, float]] = None):
        self.weights = weights or DEFAULT_STEMNESS_WEIGHTS

    def predict(self, expression_df: pd.DataFrame) -> pd.Series:
        """
        Belirtilen gen ekspresyon matrisindeki her bir hücre/örnek için
        kök hücrelilik indeksi (mRNAsI) değerini [0, 1] aralığında hesaplar.

        Args:
            expression_df: [n_samples x n_genes] ekspresyon tablosu.

        Returns:
            pd.Series: Her örnek için normalize edilmiş kök hücrelilik indeksi.
        """
        common_genes = [g for g in self.weights if g in expression_df.columns]
        if not common_genes:
            raise ValueError("Ekspresyon matrisinde kök hücrelilik modeliyle eşleşen gen bulunamadı.")

        weight_vec = np.array([self.weights[g] for g in common_genes])
        sub_matrix = expression_df[common_genes].values

        # Her hücre için ağırlıklı Spearman rank korelasyonu veya skaler çarpım
        raw_scores = []
        for i in range(sub_matrix.shape[0]):
            sample_expr = sub_matrix[i, :]
            # Sıralama korelasyonu
            corr, _ = spearmanr(sample_expr, weight_vec)
            if np.isnan(corr):
                # Sıfır varyans durumunda skaler çarpım
                corr = np.dot(sample_expr, weight_vec) / (np.linalg.norm(sample_expr) * np.linalg.norm(weight_vec) + 1e-8)
            raw_scores.append(corr)

        raw_scores = np.array(raw_scores)

        # Min-Max normalizasyonu ile [0, 1] ölçeğine taşıma (Malta et al. metodu)
        min_val = np.min(raw_scores)
        max_val = np.max(raw_scores)
        if max_val > min_val:
            normalized_scores = (raw_scores - min_val) / (max_val - min_val)
        else:
            normalized_scores = np.ones_like(raw_scores) * 0.5

        return pd.Series(normalized_scores, index=expression_df.index, name="stemness_index")


def calculate_stemness_index(expression_df: pd.DataFrame) -> pd.Series:
    """Tek fonksiyonla kök hücrelilik indeksini döndüren pratik arayüz."""
    model = StemnessModel()
    return model.predict(expression_df)
