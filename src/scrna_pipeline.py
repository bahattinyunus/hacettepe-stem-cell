"""
Tek Hücre RNA Dizileme (scRNA-seq) Analiz Hattı ve Kök Hücre Veri Üreteci.

Tek hücre seviyesinde kalite kontrolü, normalizasyon, değişken gen seçimi,
boyutsallık indirgeme (PCA / UMAP projeksiyonu), marker gen analizi ve
pseudotime (hücre kaderi zamanı) modellemesi sağlar.
"""

from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


# Temel kök hücre, MSC, podosit ve CSC belirteç genleri
CURATED_MARKER_GENES = {
    "Pluripotensi": ["POU5F1", "SOX2", "NANOG", "KLF4", "MYC", "ZFP42"],
    "Somatik/Fibroblast": ["COL1A1", "COL1A2", "ACTA2", "FAP"],
    "MSC_Pozitif": ["NT5E", "THY1", "ENG", "VCAM1"],  # CD73, CD90, CD105
    "MSC_Negatif": ["PTPRC", "CD34", "ITGAM"],         # CD45, CD34, CD11b
    "Podosit": ["NPHS1", "NPHS2", "WT1", "SYNPO"],     # Glomerüler podosit belirteçleri
    "Kanser_Kok_Hucresi": ["THY1", "PROM1", "CD44", "ALDH1A1", "LGR5", "EPCAM"],
    "Evrensel_ve_Mito": ["GAPDH", "ACTB", "MT-CO1", "MT-ND4"],
}


def generate_synthetic_stem_cell_data(
    n_cells: int = 600,
    n_genes: int = 50,
    seed: int = 42,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Hacettepe Kök Hücre Çalıştayı oturumlarındaki biyolojik senaryoları yansıtan
    sentetik tek hücre (scRNA-seq) ekspresyon matrisi ve hücre metaverisi üretir.

    Hücre Tipleri:
    1. Pluripotent ESC (Embriyonik Kök Hücre)
    2. Erken iPSC (Yamanaka faktörleriyle reprogramlanan hücreler)
    3. Geç iPSC (Pluripotensiye tam kavuşan hücreler)
    4. Somatik Fibroblast (Deri/bağ doku hücresi)
    5. Kemik İliği Mezenkimal Kök Hücresi (BM-MSC)
    6. Hasarlı Podosit (Böbrek glomerül hücresi)
    7. Kolorektal Kanser Kök Hücresi (THY1/CD90+ CSC)

    Returns:
        counts_df: [n_cells x n_genes] UMI sayı matrisi
        metadata_df: [n_cells x metadata_fields] hücre etiketleri, kütüphane büyüklüğü vb.
    """
    np.random.seed(seed)

    cell_types = [
        ("Pluripotent ESC", int(n_cells * 0.15)),
        ("Erken iPSC", int(n_cells * 0.12)),
        ("Geç iPSC", int(n_cells * 0.13)),
        ("Somatik Fibroblast", int(n_cells * 0.18)),
        ("Kemik Iligi MSC", int(n_cells * 0.16)),
        ("Hasarli Podosit", int(n_cells * 0.14)),
        ("Kolorektal CSC", n_cells - (
            int(n_cells * 0.15) + int(n_cells * 0.12) + int(n_cells * 0.13) +
            int(n_cells * 0.18) + int(n_cells * 0.16) + int(n_cells * 0.14)
        )),
    ]

    # Tüm gen listesini tekilleştirerek düzleştir
    all_curated_genes = []
    for g_list in CURATED_MARKER_GENES.values():
        for g in g_list:
            if g not in all_curated_genes:
                all_curated_genes.append(g)
    
    # Ek genler ekleyerek zenginleştir
    extra_genes = [f"GENE_{i:03d}" for i in range(max(0, n_genes - len(all_curated_genes)))]
    gene_names = all_curated_genes[:n_genes] + extra_genes

    cell_ids = []
    cell_type_labels = []
    expression_matrix = []

    for ctype, count in cell_types:
        for _ in range(count):
            cell_id = f"CELL_{len(cell_ids)+1:04d}"
            cell_ids.append(cell_id)
            cell_type_labels.append(ctype)

            # Bazal ekspresyon (Negatif binomiyal dağılım / Poisson)
            base_expr = np.random.negative_binomial(n=3, p=0.4, size=len(gene_names)).astype(float)

            # Tip-özel belirteç gen up-regülasyonları
            gene_dict = {g: base_expr[i] for i, g in enumerate(gene_names)}

            if ctype in ("Pluripotent ESC", "Geç iPSC"):
                for g in ["POU5F1", "SOX2", "NANOG", "KLF4", "ZFP42"]:
                    if g in gene_dict:
                        gene_dict[g] += np.random.gamma(shape=12, scale=4)
            elif ctype == "Erken iPSC":
                for g in ["KLF4", "MYC"]:
                    if g in gene_dict:
                        gene_dict[g] += np.random.gamma(shape=8, scale=3)
                for g in ["POU5F1", "NANOG"]:
                    if g in gene_dict:
                        gene_dict[g] += np.random.gamma(shape=3, scale=2)
            elif ctype == "Somatik Fibroblast":
                for g in ["COL1A1", "COL1A2", "ACTA2", "FAP"]:
                    if g in gene_dict:
                        gene_dict[g] += np.random.gamma(shape=14, scale=3)
            elif ctype == "Kemik Iligi MSC":
                for g in ["NT5E", "THY1", "ENG", "VCAM1"]:
                    if g in gene_dict:
                        gene_dict[g] += np.random.gamma(shape=15, scale=4)
                for g in ["PTPRC", "CD34"]:
                    if g in gene_dict:
                        gene_dict[g] = np.random.poisson(lam=0.2)
            elif ctype == "Hasarli Podosit":
                for g in ["NPHS1", "NPHS2", "WT1", "SYNPO"]:
                    if g in gene_dict:
                        # Hasarlı podositte nphs1 ekspresyonu düşüktür
                        gene_dict[g] += np.random.gamma(shape=4, scale=2)
            elif ctype == "Kolorektal CSC":
                for g in ["THY1", "PROM1", "CD44", "ALDH1A1", "LGR5", "EPCAM"]:
                    if g in gene_dict:
                        gene_dict[g] += np.random.gamma(shape=16, scale=4)

            # Evrensel genler
            for g in ["GAPDH", "ACTB"]:
                if g in gene_dict:
                    gene_dict[g] += np.random.gamma(shape=20, scale=3)
            # Mitokondriyal genler (QC için)
            for g in ["MT-CO1", "MT-ND4"]:
                if g in gene_dict:
                    gene_dict[g] += np.random.gamma(shape=4, scale=2)

            expression_matrix.append([gene_dict[g] for g in gene_names])

    counts_df = pd.DataFrame(expression_matrix, index=cell_ids, columns=gene_names)
    
    # Metadata tablosu
    total_counts = counts_df.sum(axis=1)
    detected_genes = (counts_df > 0).sum(axis=1)
    mito_genes = [g for g in gene_names if g.startswith("MT-")]
    mito_pct = (counts_df[mito_genes].sum(axis=1) / total_counts) * 100 if mito_genes else 0.0

    metadata_df = pd.DataFrame({
        "cell_id": cell_ids,
        "cell_type": cell_type_labels,
        "total_counts": total_counts,
        "detected_genes": detected_genes,
        "mito_percent": mito_pct,
    }, index=cell_ids)

    return counts_df, metadata_df


class ScRNAPipeline:
    """
    scRNA-seq veri ön işleme, boyutsallık indirgeme,
    kümeleme ve biyobelirteç tespit hattı.
    """

    def __init__(self, target_sum: float = 1e4):
        self.target_sum = target_sum
        self.normalized_df: Optional[pd.DataFrame] = None
        self.log1p_df: Optional[pd.DataFrame] = None
        self.hvg_genes: Optional[List[str]] = None
        self.pca_model: Optional[PCA] = None
        self.pca_coordinates: Optional[pd.DataFrame] = None
        self.pseudotime_scores: Optional[pd.Series] = None

    def preprocess(
        self,
        counts_df: pd.DataFrame,
        min_counts: int = 50,
        min_genes: int = 10,
    ) -> pd.DataFrame:
        """
        Kalite kontrol filtresi, kütüphane boyutu ölçeklemesi ve log1p dönüşümü.
        """
        # Filtreleme
        cell_mask = (counts_df.sum(axis=1) >= min_counts) & ((counts_df > 0).sum(axis=1) >= min_genes)
        filtered = counts_df.loc[cell_mask]

        # Normalizasyon: (Counts / Total_per_cell) * target_sum
        total_per_cell = filtered.sum(axis=1).values[:, np.newaxis]
        total_per_cell[total_per_cell == 0] = 1.0
        normalized = (filtered.values / total_per_cell) * self.target_sum
        self.normalized_df = pd.DataFrame(normalized, index=filtered.index, columns=filtered.columns)

        # Log dönüşümü
        self.log1p_df = np.log1p(self.normalized_df)
        return self.log1p_df

    def find_highly_variable_genes(self, top_n: int = 30) -> List[str]:
        """
        Genlerin ortalama ve varyans dağılımına göre en değişken genleri (HVG) seçer.
        """
        if self.log1p_df is None:
            raise ValueError("Önce preprocess() metodunu çalıştırmalısınız.")

        means = self.log1p_df.mean(axis=0)
        variances = self.log1p_df.var(axis=0)
        # Poisson benzeri beklenen varyans düzeltmesi
        dispersion = variances / (means + 1e-6)

        sorted_genes = dispersion.sort_values(ascending=False).index.tolist()
        self.hvg_genes = sorted_genes[:top_n]
        return self.hvg_genes

    def run_pca(self, n_components: int = 5) -> pd.DataFrame:
        """
        HVG genleri üzerinde PCA çalıştırır ve hücre koordinatlarını üretir.
        """
        if self.log1p_df is None:
            raise ValueError("Önce preprocess() metodunu çalıştırmalısınız.")

        genes_to_use = self.hvg_genes if self.hvg_genes is not None else self.log1p_df.columns.tolist()
        data = self.log1p_df[genes_to_use].values

        scaler = StandardScaler()
        scaled_data = scaler.fit_transform(data)

        self.pca_model = PCA(n_components=n_components, random_state=42)
        coords = self.pca_model.fit_transform(scaled_data)

        cols = [f"PC_{i+1}" for i in range(n_components)]
        self.pca_coordinates = pd.DataFrame(coords, index=self.log1p_df.index, columns=cols)
        return self.pca_coordinates

    def cluster_cells(self, n_clusters: int = 5) -> pd.Series:
        """K-means algoritması ile hücreleri fenotipik kümelere ayırır."""
        if self.pca_coordinates is None:
            self.run_pca()
        kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
        clusters = kmeans.fit_predict(self.pca_coordinates.values)
        return pd.Series(clusters, index=self.pca_coordinates.index, name="cluster")

    def compute_pseudotime(self, root_cell_id: Optional[str] = None) -> pd.Series:
        """
        Kök hücreden farklılaşma yönüne doğru pseudotime (hücre gelişim zamanı) skoru hesaplar.
        Öklid mesafesi ve PC1/PC2 izdüşümü kullanılır.
        """
        if self.pca_coordinates is None:
            self.run_pca()

        coords = self.pca_coordinates[["PC_1", "PC_2"]].values

        if root_cell_id and root_cell_id in self.pca_coordinates.index:
            root_idx = self.pca_coordinates.index.get_loc(root_cell_id)
            root_point = coords[root_idx]
        else:
            # Otomatik: Minimum PC1 değerine sahip hücre (Pluripotent kök hücresi kabul edilir)
            root_point = coords[np.argmin(coords[:, 0])]

        distances = np.linalg.norm(coords - root_point, axis=1)
        # 0 ile 1 arasına normalize et
        pt_normalized = (distances - distances.min()) / (distances.max() - distances.min() + 1e-8)
        self.pseudotime_scores = pd.Series(pt_normalized, index=self.pca_coordinates.index, name="pseudotime")
        return self.pseudotime_scores

    def differential_expression_gene(
        self,
        gene: str,
        group_a_cells: List[str],
        group_b_cells: List[str],
    ) -> Dict[str, float]:
        """İki hücre grubu arasında belirli bir genin kat artışını ve p-değerini hesaplar."""
        if self.log1p_df is None:
            raise ValueError("Önce preprocess() metodunu çalıştırmalısınız.")

        vals_a = self.log1p_df.loc[group_a_cells, gene].values
        vals_b = self.log1p_df.loc[group_b_cells, gene].values

        mean_a = float(np.mean(vals_a))
        mean_b = float(np.mean(vals_b))
        log2_fc = float(mean_a - mean_b)  # log1p ölçeğinde yaklaşık logFC

        ttest_res = stats.ttest_ind(vals_a, vals_b, equal_var=False)
        p_val = float(ttest_res.pvalue) if not np.isnan(ttest_res.pvalue) else 1.0

        return {
            "gene": gene,
            "mean_group_a": mean_a,
            "mean_group_b": mean_b,
            "log2_fold_change": log2_fc,
            "p_value": p_val,
        }
