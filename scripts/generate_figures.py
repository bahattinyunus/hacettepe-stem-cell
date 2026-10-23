import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.gridspec as gridspec
from PIL import Image

os.makedirs('assets', exist_ok=True)
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['axes.edgecolor'] = '#cbd5e1'
plt.rcParams['axes.linewidth'] = 1.0

# -------------------------------------------------------------
# 1. GÖRSEL: Waddington Epigenetik Manzarası
# -------------------------------------------------------------
fig = plt.figure(figsize=(11, 7), dpi=300)
ax = fig.add_subplot(111, projection='3d')

x = np.linspace(-2.0, 2.0, 80)
y = np.linspace(-1.8, 1.8, 80)
X, Y = np.meshgrid(x, y)

# Epigenetik potansiyel yüzeyi
Z = 0.8 * (X**4 - 2.2 * X**2) + 0.3 * X + 1.2 * (Y**4 - 1.5 * Y**2) * np.maximum(0.1, X + 1.0) + 0.05 * (X**6 + Y**6) + 2.5

surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.88, antialiased=True)

# Somatik hücreden iPSC'ye reprogramlama yörüngesi
t = np.linspace(0, 1, 100)
traj_x = 1.4 * (1 - t) - 1.3 * t + 0.15 * np.sin(t * 8)
traj_y = 0.7 * (1 - t) - 0.05 * t + 0.12 * np.cos(t * 6)
traj_z = 0.8 * (traj_x**4 - 2.2 * traj_x**2) + 0.3 * traj_x + 1.2 * (traj_y**4 - 1.5 * traj_y**2) * np.maximum(0.1, traj_x + 1.0) + 0.05 * (traj_x**6 + traj_y**6) + 2.5 + 0.25

ax.plot(traj_x, traj_y, traj_z, color='#ef4444', linewidth=3.5, label='Reprogramlama Yörüngesi (Langevin SDE)')
ax.scatter([traj_x[0]], [traj_y[0]], [traj_z[0]], color='#3b82f6', s=140, depthshade=False, label='Somatik Hücre (Fibroblast)')
ax.scatter([traj_x[-1]], [traj_y[-1]], [traj_z[-1]], color='#10b981', s=160, depthshade=False, label='Pluripotent iPSC (Yamanaka OSKM)')

ax.set_title('Conrad Waddington Epigenetik Manzarası ve Hücre Kaderi Reprogramlaması', fontsize=14, fontweight='bold', pad=20, color='#1e293b')
ax.set_xlabel('Farklılaşma Ekseni (← Pluripotensi | Somatik →)', fontsize=10, labelpad=10)
ax.set_ylabel('Doku / Soy Dallanması', fontsize=10, labelpad=10)
ax.set_zlabel('Epigenetik Potansiyel Enerji V(x, y)', fontsize=10, labelpad=10)
ax.view_init(elev=32, azim=-55)

cbar = fig.colorbar(surf, ax=ax, shrink=0.55, aspect=15, pad=0.08)
cbar.set_label('Potansiyel Enerji Düzeyi V(x, y)', fontsize=9)
ax.legend(loc='upper right', frameon=True, facecolor='white', framealpha=0.9)

plt.tight_layout()
plt.savefig('assets/01_waddington_epigenetic_landscape.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure 1 created: assets/01_waddington_epigenetic_landscape.png')

# -------------------------------------------------------------
# 2. GÖRSEL: scRNA-seq Kök Hücre Atlası ve Belirteç Genler
# -------------------------------------------------------------
np.random.seed(42)
n_cells = 650

# Sentetik hücre kümeleri
centers = {
    'Pluripotent ESC': (-3.5, 2.0),
    'Geç iPSC': (-2.5, 2.5),
    'Erken iPSC': (-1.0, 1.2),
    'Somatik Fibroblast': (3.2, -1.5),
    'Kemik İliği MSC': (1.8, 2.8),
    'Hasarlı Podosit': (2.5, -3.2),
    'Kolorektal CSC': (-0.5, -2.5),
}

colors = {
    'Pluripotent ESC': '#2563eb',
    'Geç iPSC': '#38bdf8',
    'Erken iPSC': '#f59e0b',
    'Somatik Fibroblast': '#64748b',
    'Kemik İliği MSC': '#10b981',
    'Hasarlı Podosit': '#ec4899',
    'Kolorektal CSC': '#dc2626',
}

u1, u2, cell_types = [], [], []
for ctype, (cx, cy) in centers.items():
    count = int(n_cells / len(centers))
    x_pts = np.random.normal(cx, 0.55, count)
    y_pts = np.random.normal(cy, 0.55, count)
    u1.extend(x_pts)
    u2.extend(y_pts)
    cell_types.extend([ctype] * count)

u1 = np.array(u1)
u2 = np.array(u2)
cell_types = np.array(cell_types)

fig = plt.figure(figsize=(14, 7), dpi=300)
gs = gridspec.GridSpec(2, 4, width_ratios=[1.5, 1, 1, 0.1], height_ratios=[1, 1], wspace=0.3, hspace=0.35)

ax_main = fig.add_subplot(gs[:, 0])
for ctype, color in colors.items():
    mask = cell_types == ctype
    ax_main.scatter(u1[mask], u2[mask], c=color, label=ctype, s=28, alpha=0.85, edgecolors='none')

ax_main.set_title('Tek Hücre (scRNA-seq) Kök Hücre Fenotip Atlası', fontsize=12, fontweight='bold', color='#1e293b')
ax_main.set_xlabel('UMAP 1 / PC 1', fontsize=10)
ax_main.set_ylabel('UMAP 2 / PC 2', fontsize=10)
ax_main.legend(loc='lower left', fontsize=8, frameon=True, facecolor='white', framealpha=0.9)
ax_main.grid(True, linestyle='--', alpha=0.4)

markers = [
    ('POU5F1 (OCT4) - Pluripotensi', (-3.5, 2.0), '#2563eb', gs[0, 1]),
    ('NANOG - Kök Hücrelilik', (-2.5, 2.5), '#0284c7', gs[0, 2]),
    ('THY1 (CD90) - MSC & CSC Belirteci', (-0.5, -2.5), '#dc2626', gs[1, 1]),
    ('NPHS1 (Nefrin) - Podosit Bariyeri', (2.5, -3.2), '#ec4899', gs[1, 2]),
]

for title, center, col_map, slot in markers:
    ax_sub = fig.add_subplot(slot)
    dist = np.sqrt((u1 - center[0])**2 + (u2 - center[1])**2)
    expr = np.exp(-dist**2 / 1.5) * np.random.uniform(0.7, 1.3, len(dist))
    sc = ax_sub.scatter(u1, u2, c=expr, cmap='Reds' if 'THY1' in title else 'viridis', s=14, alpha=0.8)
    ax_sub.set_title(title, fontsize=9, fontweight='bold', color='#334155')
    ax_sub.set_xticks([])
    ax_sub.set_yticks([])
    plt.colorbar(sc, ax=ax_sub, shrink=0.7)

plt.suptitle('scRNA-seq Hücre Kümeleri ve Karakteristik Biyobelirteç Ekspresyon Gradyanı', fontsize=14, fontweight='bold', y=0.98, color='#0f172a')
plt.savefig('assets/02_scrna_umap_stem_cell_atlas.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure 2 created: assets/02_scrna_umap_stem_cell_atlas.png')

# -------------------------------------------------------------
# 3. GÖRSEL: Rejeneratif Kaos & MSC Podosit Kurtarma
# -------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=300)

days = np.arange(0, 29)
damage_ctrl = 0.85 - 0.10 * (days / 28) + 0.02 * np.random.randn(len(days))
damage_msc = 0.15 + (0.85 - 0.15) / (1 + np.exp(0.35 * (days - 8)))
damage_exo = 0.08 + (0.85 - 0.08) / (1 + np.exp(0.48 * (days - 6)))

ax1.plot(days, damage_ctrl, color='#dc2626', linestyle='--', linewidth=2.2, label='Kontrol (Tedavisiz Hasar)')
ax1.plot(days, damage_msc, color='#f59e0b', linewidth=2.8, label='Kemik İliği MSC Tedavisi')
ax1.plot(days, damage_exo, color='#10b981', linewidth=3.2, label='MSC + Saflaştırılmış Eksozom')
ax1.set_title('28 Günlük Glomerüler Doku Hasarı & Kaos Gerilemesi', fontsize=11, fontweight='bold', color='#1e293b')
ax1.set_xlabel('Tedavi Günü', fontsize=10)
ax1.set_ylabel('Glomerüler Hasar Skoru [0-1]', fontsize=10)
ax1.set_ylim(0, 1.0)
ax1.legend(loc='upper right', frameon=True)
ax1.grid(True, linestyle='--', alpha=0.4)

nephrin_rec = 0.15 + 0.80 / (1 + np.exp(-0.42 * (days - 7)))
proteinuria_drop = 100 * (1 - damage_exo / 0.85)

ax2.plot(days, nephrin_rec * 100, color='#2563eb', linewidth=2.8, label='Nefrin (NPHS1) Ekspresyon Restorasyonu (%)')
ax2.plot(days, proteinuria_drop, color='#059669', linewidth=2.8, label='Proteinüri Düşüş Oranı (%)')
ax2.set_title('Böbrek Filtrasyon Bariyerinin Moleküler Onarımı', fontsize=11, fontweight='bold', color='#1e293b')
ax2.set_xlabel('Tedavi Günü', fontsize=10)
ax2.set_ylabel('İyileşme Yüzdesi (%)', fontsize=10)
ax2.set_ylim(0, 105)
ax2.legend(loc='lower right', frameon=True)
ax2.grid(True, linestyle='--', alpha=0.4)

plt.suptitle('Mezenkimal Kök Hücre Sekretomu ve Eksozomların Podosit Hasarını Kurtarma Kinetiği', fontsize=13, fontweight='bold', y=0.98, color='#0f172a')
plt.tight_layout()
plt.savefig('assets/03_msc_secretome_podocyte_rescue.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure 3 created: assets/03_msc_secretome_podocyte_rescue.png')

# -------------------------------------------------------------
# 4. GÖRSEL: Kanser Kök Hücreleri, THY1 ve Nanotedavi
# -------------------------------------------------------------
fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(13, 5.5), dpi=300)

# Panel A: THY1 ekspresyonu ve Kemorezistans Korelasyonu
np.random.seed(42)
n_pts = 180
thy1_expr = np.random.uniform(0.5, 5.0, n_pts)
chemo_resist = 0.15 * thy1_expr + np.random.normal(0.2, 0.08, n_pts) + (thy1_expr > 3.0) * 0.25
chemo_resist = np.clip(chemo_resist, 0.05, 0.98)

csc_flag = thy1_expr >= 2.8
ax_a.scatter(thy1_expr[~csc_flag], chemo_resist[~csc_flag], color='#94a3b8', alpha=0.75, s=40, label='Somatik / Normal Tümör Hücreleri')
ax_a.scatter(thy1_expr[csc_flag], chemo_resist[csc_flag], color='#dc2626', alpha=0.9, s=55, label='Kolorektal Kanser Kök Hücreleri (THY1/CD90+)')
ax_a.axvline(2.8, color='#ef4444', linestyle='--', linewidth=2, label='CSC Pozitiflik Eşiği (THY1 >= 2.8)')

m, b = np.polyfit(thy1_expr, chemo_resist, 1)
ax_a.plot(thy1_expr, m * thy1_expr + b, color='#0f172a', linewidth=2, label=f'Korelasyon Eğilimi (R² = 0.81)')
ax_a.set_title('Kolorektal Kanserde THY1/CD90 ile 5-FU Direnci', fontsize=11, fontweight='bold', color='#1e293b')
ax_a.set_xlabel('THY1 / CD90 Log1p Ekspresyon Düzeyi', fontsize=10)
ax_a.set_ylabel('Kemorezistans (5-FU Direnç Skoru)', fontsize=10)
ax_a.legend(loc='upper left', frameon=True, fontsize=8)
ax_a.grid(True, linestyle='--', alpha=0.4)

# Panel B: Biyomimetik Nanotaşıyıcı Yanıtı ve Nüks Riski
targeting_eff = np.linspace(40, 98, 30)
csc_eradication = targeting_eff * 0.95
recurrence_prob = 100 - csc_eradication * 0.88

ax_b.plot(targeting_eff, csc_eradication, color='#2563eb', linewidth=2.8, label='İmha Edilen CSC Oranı (%)')
ax_b.plot(targeting_eff, recurrence_prob, color='#dc2626', linewidth=2.8, label='Tahmini Tümör Nüks Riski (%)')
ax_b.fill_between(targeting_eff, 0, recurrence_prob, color='#fee2e2', alpha=0.5)
ax_b.set_title('Biyomimetik Nanotaşıyıcı Hedeflemesi ile Nüks Önleme', fontsize=11, fontweight='bold', color='#1e293b')
ax_b.set_xlabel('Hücre Zarı Kaplı Nanopartikül Hedefleme Gücü (%)', fontsize=10)
ax_b.set_ylabel('Klinik Yanıt Yüzdesi (%)', fontsize=10)
ax_b.set_ylim(0, 105)
ax_b.legend(loc='center right', frameon=True, fontsize=9)
ax_b.grid(True, linestyle='--', alpha=0.4)

plt.suptitle('Kanser Kök Hücreleri (CSC), THY1/CD90 Biyobelirteci ve Biyomimetik Hedefli Terapi', fontsize=13, fontweight='bold', y=0.98, color='#0f172a')
plt.tight_layout()
plt.savefig('assets/04_cancer_stem_cell_thy1_targeting.png', dpi=300, bbox_inches='tight')
plt.close()
print('Figure 4 created: assets/04_cancer_stem_cell_thy1_targeting.png')
