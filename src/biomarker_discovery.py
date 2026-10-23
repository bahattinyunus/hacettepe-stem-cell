"""
Multi-Omik Entegrasyon ve Biyobelirteç Keşfi (Biomarker Discovery) Modülü.

"Binlerce Sinyal, Kaç Biyobelirteç?"
Transkriptomik ve proteomik sinyaller arasından istatistiksel diferansiyel analiz
(Volcano Plot, Benjamini-Hochberg FDR) ve Makine Öğrenmesi (Random Forest, Lasso)
kullanarak minimal ve yüksek doğruluklu biyobelirteç panelleri seçer.
"""

from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_curve, auc, classification_report, accuracy_score


class DifferentialExpression:
    """
    İki grup arasında çoklu hipotez testi ve kat artışı analizi yapar.
    """

    @staticmethod
    def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
        """Benjamini-Hochberg FDR (False Discovery Rate) düzeltmesi."""
        n = len(p_values)
        sorted_indices = np.argsort(p_values)
        sorted_p = p_values[sorted_indices]

        q_values = np.zeros(n)
        cummin = 1.0

        for i in range(n - 1, -1, -1):
            rank = i + 1
            q = sorted_p[i] * n / rank
            cummin = min(cummin, q)
            q_values[i] = min(cummin, 1.0)

        # Orijinal sıraya geri döndür
        rev_indices = np.empty(n, dtype=int)
        rev_indices[sorted_indices] = np.arange(n)
        return q_values[rev_indices]

    @classmethod
    def analyze(
        cls,
        expression_df: pd.DataFrame,
        group_a_samples: List[str],
        group_b_samples: List[str],
        log2_fc_threshold: float = 1.0,
        fdr_threshold: float = 0.05,
    ) -> pd.DataFrame:
        """
        Her gen için Log2 Fold Change, P-değeri, Q-değeri ve regülasyon durumunu belirler.
        """
        genes = expression_df.columns
        res = []

        mat_a = expression_df.loc[group_a_samples].values
        mat_b = expression_df.loc[group_b_samples].values

        p_vals = []
        log2_fcs = []
        mean_as = []
        mean_bs = []

        for j in range(len(genes)):
            vals_a = mat_a[:, j]
            vals_b = mat_b[:, j]

            mean_a = float(np.mean(vals_a))
            mean_b = float(np.mean(vals_b))

            # Log2 Fold Change (küçük epsilon ile sıfıra bölmeyi engelle)
            fc = (mean_a + 1e-4) / (mean_b + 1e-4)
            log2_fc = float(np.log2(fc))

            ttest = stats.ttest_ind(vals_a, vals_b, equal_var=False)
            pval = float(ttest.pvalue) if not np.isnan(ttest.pvalue) else 1.0

            mean_as.append(mean_a)
            mean_bs.append(mean_b)
            log2_fcs.append(log2_fc)
            p_vals.append(pval)

        p_vals = np.array(p_vals)
        q_vals = cls.benjamini_hochberg(p_vals)
        neg_log10_p = -np.log10(np.maximum(p_vals, 1e-300))

        # Regülasyon sınıflandırması
        status = []
        for lfc, q in zip(log2_fcs, q_vals):
            if q < fdr_threshold and lfc >= log2_fc_threshold:
                status.append("Up-regulated")
            elif q < fdr_threshold and lfc <= -log2_fc_threshold:
                status.append("Down-regulated")
            else:
                status.append("Not Significant")

        df = pd.DataFrame({
            "gene": genes,
            "mean_group_a": mean_as,
            "mean_group_b": mean_bs,
            "log2_fold_change": log2_fcs,
            "p_value": p_vals,
            "fdr_q_value": q_vals,
            "-log10_p_value": neg_log10_p,
            "status": status,
        })
        return df.sort_values(by="fdr_q_value")


class BiomarkerClassifier:
    """
    Makine öğrenmesi ile biyobelirteç seçimi ve sınıflandırma modeli.
    """

    def __init__(self, n_estimators: int = 100, random_state: int = 42):
        self.rf = RandomForestClassifier(n_estimators=n_estimators, random_state=random_state)
        self.feature_importances_: Optional[pd.Series] = None
        self.selected_biomarkers: Optional[List[str]] = None

    def fit_and_rank(self, X: pd.DataFrame, y: pd.Series, top_k: int = 10) -> pd.Series:
        """
        Öznitelik önem derecelerini hesaplar ve en ayırt edici 'top_k' biyobelirteci seçer.
        """
        self.rf.fit(X, y)
        importances = pd.Series(self.rf.feature_importances_, index=X.columns)
        self.feature_importances_ = importances.sort_values(ascending=False)
        self.selected_biomarkers = self.feature_importances_.head(top_k).index.tolist()
        return self.feature_importances_

    def evaluate_cv(
        self,
        X: pd.DataFrame,
        y: pd.Series,
        n_splits: int = 5,
    ) -> Dict[str, Union[float, np.ndarray, Dict]]:
        """
        Stratified K-Fold Çapraz Doğrulama ile ROC-AUC ve doğruluk hesaplar.
        """
        cv = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
        
        # Olasılık tahminleri
        y_probs = cross_val_predict(self.rf, X, y, cv=cv, method="predict_proba")[:, 1]
        y_preds = cross_val_predict(self.rf, X, y, cv=cv)

        fpr, tpr, _ = roc_curve(y, y_probs)
        roc_auc = float(auc(fpr, tpr))
        accuracy = float(accuracy_score(y, y_preds))
        report = classification_report(y, y_preds, output_dict=True)

        return {
            "roc_auc": roc_auc,
            "accuracy": accuracy,
            "fpr": fpr,
            "tpr": tpr,
            "y_probs": y_probs,
            "report": report,
        }
