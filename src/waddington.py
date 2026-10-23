"""
Waddington Epigenetik Manzarası ve Yamanaka Yeniden Programlama Simülatörü.

Conrad Waddington'ın 1957 yılında öne sürdüğü epigenetik manzara modeli,
hücre kaderi ve farklılaşma vadilerini potansiyel enerji yüzeyi üzerinde modeller.
Yamanaka faktörleri (Oct4, Sox2, Klf4, c-Myc) bu bariyeri düşürerek somatik hücrelerin
pluripotent kök hücre vadisine geri dönmesini (reprogramming) sağlar.
"""

from typing import Dict, List, Optional, Tuple
import numpy as np


class WaddingtonLandscape:
    """
    Epigenetik enerji potansiyeli fonksiyonunu hesaplayan ve
    hücresel dinamikleri simüle eden sınıf.
    """

    def __init__(self, barrier_height: float = 2.5, asymmetry: float = 0.3):
        """
        Args:
            barrier_height: Pluripotensi ve somatik durum arasındaki potansiyel bariyer yüksekliği.
            asymmetry: Farklılaşma yönündeki asimetri eğimi (Waddington meyli).
        """
        self.barrier_height = barrier_height
        self.asymmetry = asymmetry

    def potential(self, x: np.ndarray, y: np.ndarray, yamanaka_factor_strength: float = 0.0) -> np.ndarray:
        """
        Belirtilen koordinatlarda Waddington potansiyel enerjisi V(x, y).

        x ekseni: Farklılaşma koordinatı (-1.5: Pluripotent ESC/iPSC, +1.5: Somatik hücre)
        y ekseni: Lineage / Soy kollarındaki dallanma (-1.0: Mezenkimal/Stromal, +1.0: Epitelyal/Podosit)
        yamanaka_factor_strength: Reprogramlama faktörlerinin konsantrasyonu [0.0 - 1.0]

        Matematiksel formülasyon:
        V(x, y) = a*(x^4 - 2*x^2) + b*y^2*(x + 1.2) + c*x - yamanaka * driving_force
        """
        a = 0.8 * (1.0 - 0.7 * yamanaka_factor_strength)
        barrier = self.barrier_height * (1.0 - 0.85 * yamanaka_factor_strength)
        x_pot = a * (x**4 - 2.2 * (x**2)) + (self.asymmetry - 1.2 * yamanaka_factor_strength) * x

        # Dallanma vadileri: x arttıkça (farklılaştıkça) y eksenindeki iki tepe arası açılır
        branching = 1.2 * (y**4 - 1.5 * (y**2)) * np.maximum(0.1, x + 1.0)
        
        # Sınır cezası
        boundary = 0.05 * (x**6 + y**6)

        return x_pot + branching + boundary + barrier

    def gradient(self, x: float, y: float, yamanaka_factor_strength: float = 0.0, eps: float = 1e-4) -> Tuple[float, float]:
        """Potansiyelin sayısal gradyanı (-dV/dx, -dV/dy)."""
        v_xp = self.potential(np.array([x + eps]), np.array([y]), yamanaka_factor_strength)[0]
        v_xm = self.potential(np.array([x - eps]), np.array([y]), yamanaka_factor_strength)[0]
        grad_x = (v_xp - v_xm) / (2 * eps)

        v_yp = self.potential(np.array([x]), np.array([y + eps]), yamanaka_factor_strength)[0]
        v_ym = self.potential(np.array([x]), np.array([y - eps]), yamanaka_factor_strength)[0]
        grad_y = (v_yp - v_ym) / (2 * eps)

        return float(grad_x), float(grad_y)

    def generate_surface_grid(
        self,
        grid_points: int = 50,
        x_range: Tuple[float, float] = (-2.0, 2.0),
        y_range: Tuple[float, float] = (-1.8, 1.8),
        yamanaka_factor_strength: float = 0.0,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Plotly ve 3D çizimler için X, Y, Z yüzey matrislerini üretir."""
        x = np.linspace(x_range[0], x_range[1], grid_points)
        y = np.linspace(y_range[0], y_range[1], grid_points)
        xx, yy = np.meshgrid(x, y)
        zz = self.potential(xx, yy, yamanaka_factor_strength=yamanaka_factor_strength)
        return xx, yy, zz


def simulate_reprogramming(
    start_pos: Tuple[float, float] = (1.4, 0.8),
    yamanaka_strength: float = 0.85,
    steps: int = 250,
    dt: float = 0.04,
    diffusion_noise: float = 0.12,
    seed: Optional[int] = 42,
) -> Dict[str, np.ndarray]:
    """
    Langevin stokastik diferansiyel denklemi ile somatik hücreden iPSC'ye geçiş yörüngesini simüle eder:
    dX_t = -grad(V(X_t)) * dt + sigma * dW_t

    Args:
        start_pos: Somatik hücre başlangıç koordinatı (x, y)
        yamanaka_strength: Yamanaka kokteyli dozu (0.0: somatik kalır, >0.6: pluripotensiye döner)
        steps: Adım sayısı
        dt: Zaman adımı
        diffusion_noise: Biyolojik gürültü / stokastisite derecesi
        seed: Rastgelelik çekirdeği

    Returns:
        Yörünge koordinatları, potansiyel enerjisi ve hücre durumu etiketleri.
    """
    if seed is not None:
        np.random.seed(seed)

    landscape = WaddingtonLandscape()
    trajectory = np.zeros((steps, 2))
    energy = np.zeros(steps)
    states = []

    curr_x, curr_y = start_pos

    for step in range(steps):
        # Zamanla yamanaka faktörlerinin ekspresyon dinamikleri
        eff_yamanaka = yamanaka_strength * min(1.0, (step + 10) / (steps * 0.4))
        
        gx, gy = landscape.gradient(curr_x, curr_y, yamanaka_factor_strength=eff_yamanaka)
        
        # Langevin güncellemesi
        noise_x = np.random.normal(0, np.sqrt(dt) * diffusion_noise)
        noise_y = np.random.normal(0, np.sqrt(dt) * diffusion_noise)

        curr_x = curr_x - gx * dt + noise_x
        curr_y = curr_y - gy * dt + noise_y

        # Sınır koruması
        curr_x = np.clip(curr_x, -2.1, 2.1)
        curr_y = np.clip(curr_y, -1.9, 1.9)

        trajectory[step] = [curr_x, curr_y]
        energy[step] = landscape.potential(np.array([curr_x]), np.array([curr_y]), yamanaka_factor_strength=eff_yamanaka)[0]

        # Hücre kimliği sınıflandırması
        if curr_x < -0.8:
            state = "Pluripotent (iPSC/ESC)"
        elif curr_x < 0.2:
            state = "Ara Durum (Pre-iPSC / Epigenetik Geçiş)"
        elif curr_x < 1.0:
            state = "Kısmi Yeniden Programlanmış"
        else:
            state = "Somatik Hücre (Farklılaşmış)"
        states.append(state)

    return {
        "trajectory": trajectory,
        "energy": energy,
        "states": np.array(states),
        "steps": np.arange(steps),
        "final_state": states[-1],
        "success": bool(trajectory[-1, 0] < -0.6),
    }
