"""
🧬 In Silico Whole-Cell Digital Twin & Perturbation Engine
Hacettepe Stem Cell & Regenerative Medicine Computational Suite

Simulates genetic (CRISPR KO / Overexpression), chemical (small molecule inhibitors),
and trophic perturbations on cellular regulatory networks based on experimental findings
from the Hacettepe University Stem Cell Technologies Workshop.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd


# Default baseline expression across key phenotypic networks
BASELINE_EXPRESSION: Dict[str, float] = {
    # Pluripotency / Stemness
    "POU5F1": 0.15,
    "SOX2": 0.12,
    "NANOG": 0.10,
    # Epigenetic barriers
    "SETD2": 1.20,
    "DNMT1": 1.10,
    "MAT2A": 0.95,
    # Glomerular / Podocyte integrity
    "NPHS1": 1.50,
    "NPHS2": 1.40,
    # Inflammatory / Chemotaxis
    "MEFV": 0.80,
    "IL1B": 0.40,
    "CXCR4": 0.70,
    # Stroma / Fibrosis / CSC
    "THY1": 0.90,
    "CEBPB": 0.60,
    "PPARG": 0.80,
    "ADIPOQ": 0.75,
    # Survival / Resistance
    "BCL2": 1.00,
    "BAX": 0.50,
    "ABCG1": 0.65,
    "MMP2": 0.55,
}

# Gene interaction weight matrix (Direct activations (+) and repressions (-))
REGULATORY_MATRIX: Dict[str, Dict[str, float]] = {
    "SETD2": {"POU5F1": -0.85, "SOX2": -0.80, "NANOG": -0.75},  # H3K36me3 epigenetic barrier
    "MEFV": {"IL1B": 0.70, "CXCR4": 0.60},                      # Pyrin inflammasome axis
    "THY1": {"CXCR4": 0.65, "MMP2": 0.70, "ABCG1": 0.55},       # Stroma & CSC invasiveness
    "NPHS1": {"BAX": -0.75, "BCL2": 0.60},                      # Podocyte survival
    "CEBPB": {"PPARG": 0.40, "ADIPOQ": 0.30},                   # Early adipogenic regulator
    "PPARG": {"ADIPOQ": 0.85},                                   # Terminal adipogenesis
    "POU5F1": {"SOX2": 0.80, "NANOG": 0.85, "THY1": -0.30},     # Master pluripotency circuit
}


@dataclass
class PerturbationResult:
    """Stores the simulation outcome of an in silico cellular perturbation."""
    perturbation_name: str
    target: str
    perturbation_type: str  # 'CRISPR_KO', 'INHIBITION', 'OVEREXPRESSION', 'TROPHIC_RESCUE'
    initial_expression: Dict[str, float]
    perturbed_expression: Dict[str, float]
    fold_changes: Dict[str, float]
    stemness_shift: float
    viability_percentage: float
    phenotypic_transition: str
    reprogramming_efficiency_multiplier: float
    biomarker_alert: str


class WholeCellDigitalTwin:
    """
    In Silico Whole-Cell Digital Twin model for stem cell and regenerative medicine.
    Implements dynamic network propagation under targeted genetic and chemical perturbations.
    """

    def __init__(self, baseline_state: Optional[Dict[str, float]] = None, random_state: int = 42):
        self.rng = np.random.RandomState(random_state)
        self.state = baseline_state.copy() if baseline_state else BASELINE_EXPRESSION.copy()
        self.genes = sorted(list(self.state.keys()))

    def simulate_perturbation(
        self,
        target_gene_or_compound: str,
        perturbation_type: str = "CRISPR_KO",
        intensity: float = 0.90,
        noise_level: float = 0.05,
        steps: int = 50,
    ) -> PerturbationResult:
        """
        Simulates the dynamic propagation of a targeted perturbation across the cell network.

        Parameters
        ----------
        target_gene_or_compound : str
            Target gene (e.g. 'SETD2', 'MEFV', 'THY1') or compound ('GYY4137', 'Temsirolimus',
            'Niclosamide', 'Etoposide', 'L-methionine restriction').
        perturbation_type : str
            One of 'CRISPR_KO', 'INHIBITION', 'OVEREXPRESSION', 'TROPHIC_RESCUE'.
        intensity : float
            Strength of intervention [0.0 - 1.0].
        noise_level : float
            Biological stochastic variation.
        steps : int
            Euler integration steps.

        Returns
        -------
        PerturbationResult
            Complete multi-omics perturbation trajectory and phenotypic transition metrics.
        """
        initial_expr = self.state.copy()
        current_expr = self.state.copy()

        # Resolve compound effects to primary genetic / signaling targets
        primary_targets: Dict[str, float] = {}
        alert_msg = ""
        trans_label = ""
        reprog_mult = 1.0

        target_upper = target_gene_or_compound.upper()

        if "SETD2" in target_upper:
            # Epigenetic barrier knockout (Batman et al., 2024)
            primary_targets["SETD2"] = -intensity
            primary_targets["POU5F1"] = +0.8 * intensity
            primary_targets["NANOG"] = +0.75 * intensity
            reprog_mult = 1.0 + (14.5 * intensity)
            trans_label = "Somatik Hücre -> Yüksek Verimli Pluripotent iPSC Kolonisi"
            alert_msg = f"SETD2 inaktivasyonu H3K36me3 bariyerini kırdı. Reprogramlama verimi {reprog_mult:.1f}x katlandı!"

        elif "GYY4137" in target_upper or "H2S" in target_upper:
            # Slow H2S donor rescue in pulmonary fibrosis (Karaman et al., 2024)
            primary_targets["IL1B"] = -0.75 * intensity
            primary_targets["MMP2"] = -0.65 * intensity
            primary_targets["BCL2"] = +0.50 * intensity
            primary_targets["BAX"] = -0.60 * intensity
            trans_label = "Fibrotik Akciğer -> Fonksiyonel Restorasyon (Taurin & Redoks Denge)"
            alert_msg = "GYY4137 H2S salınımı kollajen/fibronektin sentezini durdurdu, akciğer kompliyansını restore etti."

        elif "NICLOSAMIDE" in target_upper:
            # STAT3/Wnt pathway inhibition in CD133+ CSCs (Altundag Erdogan et al., 2024)
            primary_targets["CXCR4"] = -0.80 * intensity
            primary_targets["MMP2"] = -0.75 * intensity
            primary_targets["ABCG1"] = -0.85 * intensity
            primary_targets["POU5F1"] = -0.60 * intensity
            trans_label = "Malign Kanser Kök Hücresi (CSC) -> Kemoduyarlı Kök Hücre İmhası"
            alert_msg = "Niklozamid Wnt/STAT3 aksını (AXIN2, LGR5, CYCLIN D1) bloke ederek nüks ve direnci yok etti."

        elif "TEMSIROLIMUS" in target_upper:
            # mTOR pathway inhibition in CSC metastasis model
            primary_targets["MMP2"] = -0.70 * intensity
            primary_targets["CXCR4"] = -0.65 * intensity
            primary_targets["ABCG1"] = -0.60 * intensity
            trans_label = "Kemik Metastaz İnvazyonu -> İnvazyon Blokajı ve Apoptoz"
            alert_msg = "Temsirolimus mTOR sinyalini keserek 3B biyomimetik kemik nişine invazyonu durdurdu."

        elif "ETOPOSIDE" in target_upper:
            # Topoisomerase II inhibitor genotoxic stress (Muratoglu et al., 2024)
            primary_targets["CEBPB"] = +0.85 * intensity  # early adipogenic genes upregulated
            primary_targets["PPARG"] = -0.60 * intensity  # terminal adipogenesis arrested
            primary_targets["ADIPOQ"] = -0.70 * intensity
            primary_targets["BAX"] = +0.75 * intensity
            trans_label = "Stromal Hücre -> Genotoksik Yaşlanma & Adipojenik Maturasyon Blokajı"
            alert_msg = "Etoposid erken adipogenezi (CEBPB) tetikledi ancak terminal PPARG/ADIPOQ maturasyonu kilitlendi."

        elif "MEFV" in target_upper:
            # FMF Pyrin knockdown (Akbaba et al., 2024)
            primary_targets["MEFV"] = -intensity
            primary_targets["IL1B"] = -0.80 * intensity
            primary_targets["CXCR4"] = -0.55 * intensity
            trans_label = "Otoinflamatuvar Nötrofil -> Polarizasyon & Kemotaksis Baskılanması"
            alert_msg = "Pyrin susturuldu; f-aktin öncü kenar (leading edge) polarizasyonu ve göç durduruldu."

        elif "THY1" in target_upper or "CD90" in target_upper:
            # Colorectal Stroma/CSC THY1 perturbation (Kayasli et al., 2024)
            primary_targets["THY1"] = -intensity
            primary_targets["CXCR4"] = -0.50 * intensity
            primary_targets["MMP2"] = -0.45 * intensity
            trans_label = "FAP+ CAF & Reaktif Stroma -> Normal Doku Düzeyi"
            alert_msg = "THY1 tükenmesi stromal reaktiviteyi ve vasküler mikroçevre sinyalleşmesini modüle etti."

        elif "PODOCYTE_RESCUE" in target_upper or "MSC_EXOSOME" in target_upper:
            # Podocyte rescue via induced MSC secretome (Sahan et al., 2024)
            primary_targets["NPHS1"] = +0.95 * intensity
            primary_targets["NPHS2"] = +0.90 * intensity
            primary_targets["BAX"] = -0.80 * intensity
            primary_targets["BCL2"] = +0.75 * intensity
            trans_label = "FSGS Glomerül Hasarı -> Nefrin Restorasyonu & Proteinüri Kurtarması"
            alert_msg = "İndüklenmiş BM-MSC parakrin faktörleri podosit sitoiskeletini onardı ve proteinüriyi sildi."

        elif target_upper in self.genes:
            sign = -1.0 if perturbation_type in ["CRISPR_KO", "INHIBITION"] else +1.0
            primary_targets[target_upper] = sign * intensity
            trans_label = f"Hedef Gen {target_upper} Modülasyonu"
            alert_msg = f"{target_upper} geni {perturbation_type} ile manipüle edildi."
        else:
            primary_targets[target_upper] = -intensity
            trans_label = f"Bilinmeyen Pertürbasyon ({target_upper})"
            alert_msg = f"{target_upper} üzerinde varsayılan regülasyon uygulandı."

        # Euler simulation of network propagation
        dt = 0.05
        gamma = 0.15  # natural degradation rate

        for _ in range(steps):
            derivs = {g: 0.0 for g in self.genes}

            # Direct perturbation forces
            for tgt, force in primary_targets.items():
                if tgt in derivs:
                    derivs[tgt] += force

            # Network propagation through regulatory interactions
            for src_gene, targets in REGULATORY_MATRIX.items():
                if src_gene in current_expr:
                    activity = current_expr[src_gene] - initial_expr[src_gene]
                    for tgt_gene, weight in targets.items():
                        if tgt_gene in derivs:
                            derivs[tgt_gene] += weight * activity * 0.4

            # Euler step with degradation and thermal stochasticity
            for g in self.genes:
                noise = self.rng.normal(0, noise_level)
                delta = (derivs[g] - gamma * (current_expr[g] - initial_expr[g])) * dt + noise
                current_expr[g] = float(np.clip(current_expr[g] + delta, 0.01, 10.0))

        # Compute log2 fold changes
        fold_changes = {}
        for g in self.genes:
            fc = np.log2((current_expr[g] + 1e-4) / (initial_expr[g] + 1e-4))
            fold_changes[g] = float(np.round(fc, 3))

        # Compute stemness shift
        pluri_genes = ["POU5F1", "SOX2", "NANOG"]
        init_pluri = np.mean([initial_expr[g] for g in pluri_genes])
        final_pluri = np.mean([current_expr[g] for g in pluri_genes])
        stemness_shift = float(np.round(final_pluri - init_pluri, 4))

        # Compute viability percentage based on BCL2/BAX balance
        bcl2_val = current_expr.get("BCL2", 1.0)
        bax_val = current_expr.get("BAX", 0.5)
        ratio = bcl2_val / (bax_val + 1e-4)
        viability = float(np.clip(50.0 + 35.0 * np.tanh(ratio - 1.5), 5.0, 99.5))

        return PerturbationResult(
            perturbation_name=target_gene_or_compound,
            target=target_upper,
            perturbation_type=perturbation_type,
            initial_expression=initial_expr,
            perturbed_expression=current_expr,
            fold_changes=fold_changes,
            stemness_shift=stemness_shift,
            viability_percentage=np.round(viability, 1),
            phenotypic_transition=trans_label,
            reprogramming_efficiency_multiplier=np.round(reprog_mult, 2),
            biomarker_alert=alert_msg,
        )

    def to_dataframe(self, result: PerturbationResult) -> pd.DataFrame:
        """Converts perturbation result into a structured pandas DataFrame."""
        rows = []
        for g in self.genes:
            rows.append({
                "Gene": g,
                "Initial": np.round(result.initial_expression[g], 3),
                "Perturbed": np.round(result.perturbed_expression[g], 3),
                "Log2FC": result.fold_changes[g],
                "Significant": abs(result.fold_changes[g]) >= 0.585,  # 1.5x threshold
            })
        return pd.DataFrame(rows).sort_values(by="Log2FC", ascending=False)
