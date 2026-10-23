import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch, ConnectionPatch

os.makedirs('assets', exist_ok=True)
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'

# ==============================================================================
# 1. GÖRSEL: EN TEPE İÇİN ESTETİK HERO BANNER (assets/00_hacettepe_stem_cell_banner.png)
# ==============================================================================
fig, ax = plt.subplots(figsize=(16, 5.5), dpi=300)
fig.patch.set_facecolor('#0a0f1d')
ax.set_facecolor('#0a0f1d')

# Arka plan degrade ve parıltı dalgaları
xlims = (0, 16)
ylims = (0, 5.5)
ax.set_xlim(xlims)
ax.set_ylim(ylims)
ax.axis('off')

# Arka plan dalgaları (Waddington epigenetik enerji çizgileri)
x_wave = np.linspace(0, 16, 500)
for i, offset in enumerate(np.linspace(0.5, 4.8, 12)):
    y_wave = offset + 0.35 * np.sin(x_wave * 0.8 + i * 0.4) + 0.15 * np.cos(x_wave * 1.5)
    alpha = 0.04 + 0.02 * (i % 3)
    ax.plot(x_wave, y_wave, color='#38bdf8', alpha=alpha, linewidth=1.2)

# Sol Taraf: Stilize DNA Çift Sarmalı
dna_x = np.linspace(0.8, 4.2, 55)
dna_y1 = 2.75 + 1.2 * np.sin(dna_x * 2.2)
dna_y2 = 2.75 - 1.2 * np.sin(dna_x * 2.2)

for i in range(0, len(dna_x), 3):
    ax.plot([dna_x[i], dna_x[i]], [dna_y1[i], dna_y2[i]], color='#0ea5e9', alpha=0.45, linewidth=1.5)
    ax.scatter([dna_x[i]], [dna_y1[i]], color='#38bdf8', s=45, alpha=0.9, edgecolors='white', linewidth=0.5)
    ax.scatter([dna_x[i]], [dna_y2[i]], color='#ec4899', s=45, alpha=0.9, edgecolors='white', linewidth=0.5)

ax.plot(dna_x, dna_y1, color='#0284c7', linewidth=2.5, alpha=0.85)
ax.plot(dna_x, dna_y2, color='#db2777', linewidth=2.5, alpha=0.85)

# Sağ Taraf: Hücresel Soy / Lineage ve Farklılaşma Ağı
nodes = [
    (14.8, 4.2, '#10b981', 'Pluripotent ESC'),
    (13.2, 3.4, '#38bdf8', 'Pre-iPSC'),
    (14.8, 2.6, '#f59e0b', 'BM-MSC'),
    (12.6, 2.0, '#a855f7', 'HSC Progenitor'),
    (14.9, 1.2, '#ec4899', 'Podosit'),
    (13.5, 0.8, '#64748b', 'Fibroblast'),
]
for x1, y1, c1, _ in nodes:
    for x2, y2, c2, _ in nodes:
        if (x1 != x2 or y1 != y2) and np.sqrt((x1-x2)**2 + (y1-y2)**2) < 2.2:
            ax.plot([x1, x2], [y1, y2], color='#334155', alpha=0.5, linewidth=1.0)

for x, y, col, lbl in nodes:
    ax.scatter([x], [y], color=col, s=180, alpha=0.95, edgecolors='white', linewidth=1.5)
    ax.text(x - 0.2, y + 0.22, lbl, color='#94a3b8', fontsize=7.5, fontweight='bold')

# Merkez Başlık Kutusu (Neon / Glassmorphism stili)
ax.text(8.2, 4.3, 'HACETTEPE ÜNİVERSİTESİ', color='#38bdf8', fontsize=12, fontweight='bold', ha='center')
ax.text(8.2, 3.45, 'KÖK HÜCRE & REJENERATİF TIP', color='#ffffff', fontsize=26, fontweight='heavy', ha='center')
ax.text(8.2, 2.75, 'Biyoinformatik • Epigenetik Reprogramlama • scRNA-seq • AI Çağında Akademi', color='#cbd5e1', fontsize=12, ha='center', style='italic')

# Rozetler (Pills)
pills = [
    (5.6, 1.6, 'Epigenetik Hafiza', '#1e293b', '#38bdf8'),
    (7.4, 1.6, 'Waddington SDE', '#1e293b', '#10b981'),
    (9.1, 1.6, 'Single-Cell Omics', '#1e293b', '#a855f7'),
    (10.9, 1.6, 'AI Dijital Ikiz', '#1e293b', '#f59e0b'),
]
for px, py, text, bg, border in pills:
    box = FancyBboxPatch((px - 0.75, py - 0.25), 1.5, 0.5, boxstyle="round,pad=0.1", fc=bg, ec=border, lw=1.2)
    ax.add_patch(box)
    ax.text(px, py, text, color='#f8fafc', fontsize=8.5, fontweight='bold', ha='center', va='center')

ax.text(8.2, 0.7, 'Sağlık Bilimleri Enstitüsü Kök Hücre Bilimleri Ana Bilim Dalı Çalıştayı ve Açık Araştırma Platformu', color='#64748b', fontsize=9.5, ha='center')

plt.tight_layout()
plt.savefig('assets/00_hacettepe_stem_cell_banner.png', dpi=300, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')
plt.close()
print('Banner created: assets/00_hacettepe_stem_cell_banner.png')


# ==============================================================================
# 2. GÖRSEL: VİZYON 1 - HÜCRESEL DİJİTAL İKİZ (assets/05_ai_whole_cell_digital_twin.png)
# ==============================================================================
fig, ax = plt.subplots(figsize=(14, 7), dpi=300)
fig.patch.set_facecolor('#0f172a')
ax.set_facecolor('#0f172a')
ax.set_xlim(0, 14)
ax.set_ylim(0, 7)
ax.axis('off')

# Başlık
ax.text(7.0, 6.5, 'AI-Güdümlü Bütüncül Hücresel Dijital İkiz (Whole-Cell Digital Twin Architecture)', color='#ffffff', fontsize=14, fontweight='bold', ha='center')
ax.text(7.0, 6.1, 'Hasta Hücresinden In Silico Simülasyona ve Kişiselleştirilmiş Kök Hücre Tedavisine', color='#94a3b8', fontsize=10, ha='center')

# Aşama 1: Multi-Omik Giriş Katmanı
box1 = FancyBboxPatch((0.5, 1.2), 3.2, 4.4, boxstyle="round,pad=0.2", fc='#1e293b', ec='#38bdf8', lw=1.5)
ax.add_patch(box1)
ax.text(2.1, 5.2, '1. Multi-Omik Giriş Katmanı', color='#38bdf8', fontsize=10.5, fontweight='bold', ha='center')
inputs = [
    ('Genomik (DNA Dizisi & SNV)', '#60a5fa'),
    ('Epigenomik (ATAC-seq & DNA Metilasyonu)', '#a78bfa'),
    ('Transkriptomik (scRNA-seq / Spatial)', '#34d399'),
    ('Proteomik & Metabolomik Sinyalleri', '#fbbf24'),
    ('Doku Mikroçevresi & Mekanik Sertlik', '#f87171'),
]
for i, (inp, col) in enumerate(inputs):
    ibox = FancyBboxPatch((0.8, 4.3 - i * 0.75), 2.6, 0.55, boxstyle="round,pad=0.08", fc='#0f172a', ec=col, lw=1.0)
    ax.add_patch(ibox)
    ax.text(2.1, 4.58 - i * 0.75, inp, color='#e2e8f0', fontsize=8, ha='center', va='center')

# Aşama 2: Biyolojik Temel Model (Foundation Model)
box2 = FancyBboxPatch((4.5, 1.2), 4.2, 4.4, boxstyle="round,pad=0.2", fc='#1e293b', ec='#a855f7', lw=1.5)
ax.add_patch(box2)
ax.text(6.6, 5.2, '2. Hücre Temel Modeli (scGPT/ESM-3)', color='#c084fc', fontsize=10.5, fontweight='bold', ha='center')

# Transformer Nöral Ağ Temsili
for lyr in range(4):
    for nd in range(5):
        cx = 5.2 + lyr * 0.95
        cy = 4.4 - nd * 0.65
        ax.scatter([cx], [cy], color='#818cf8', s=90, alpha=0.9, edgecolors='white', linewidth=0.8)
        if lyr < 3:
            for nd_next in range(5):
                ax.plot([cx, cx + 0.95], [cy, 4.4 - nd_next * 0.65], color='#475569', alpha=0.25, lw=0.6)

tbox = FancyBboxPatch((4.8, 1.5), 3.6, 0.65, boxstyle="round,pad=0.08", fc='#0f172a', ec='#818cf8', lw=1.0)
ax.add_patch(tbox)
ax.text(6.6, 1.82, '100+ Milyon Tek Hücre ile Eğitilmiş Evrensel Gizil Uzay (Latent Space)', color='#cbd5e1', fontsize=7.5, ha='center', va='center')

# Aşama 3: In Silico Simülatör ve Çıktı
box3 = FancyBboxPatch((9.5, 1.2), 4.0, 4.4, boxstyle="round,pad=0.2", fc='#1e293b', ec='#10b981', lw=1.5)
ax.add_patch(box3)
ax.text(11.5, 5.2, '3. Sanal Hücre Simülasyonu & Çıktı', color='#34d399', fontsize=10.5, fontweight='bold', ha='center')

sims = [
    ('In Silico Perturbasyon & Ilac Yaniti', 'Aylar yerine saniyeler icinde toksisite taramasi', '#34d399'),
    ('Hucre Kaderi & Reprogramlama Yolu', 'En kararli pluripotent rotanin tayini', '#38bdf8'),
    ('Minimal Yamanaka Kokteyli Optimizasyonu', 'Tumorlesme riskini sifirlayan hassas dozaj', '#fbbf24'),
    ('Hastaya Ozel Canli Biyobobrek / Organ', 'Rejeksiyon riski sifir otolog doku muhendisligi', '#f43f5e'),
]
for i, (h, sub, col) in enumerate(sims):
    sbox = FancyBboxPatch((9.8, 4.2 - i * 0.9), 3.4, 0.75, boxstyle="round,pad=0.08", fc='#0f172a', ec=col, lw=1.0)
    ax.add_patch(sbox)
    ax.text(11.5, 4.65 - i * 0.9, h, color=col, fontsize=8, fontweight='bold', ha='center', va='center')
    ax.text(11.5, 4.35 - i * 0.9, sub, color='#94a3b8', fontsize=7, ha='center', va='center')

# Bağlantı Okları
arrow1 = FancyArrowPatch((3.8, 3.4), (4.4, 3.4), arrowstyle='->', mutation_scale=20, color='#38bdf8', lw=2.5)
arrow2 = FancyArrowPatch((8.8, 3.4), (9.4, 3.4), arrowstyle='->', mutation_scale=20, color='#10b981', lw=2.5)
ax.add_patch(arrow1)
ax.add_patch(arrow2)

plt.tight_layout()
plt.savefig('assets/05_ai_whole_cell_digital_twin.png', dpi=300, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')
plt.close()
print('Vision 1 created: assets/05_ai_whole_cell_digital_twin.png')


# ==============================================================================
# 3. GÖRSEL: VİZYON 2 - OTONOM LABORATUVAR & GENÇLEŞME (assets/06_closed_loop_self_driving_lab.png)
# ==============================================================================
fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(14, 6.5), dpi=300)
fig.patch.set_facecolor('#0b1329')
for ax_sub in (ax_left, ax_right):
    ax_sub.set_facecolor('#0b1329')

# Sol Taraf: Otonom Kapalı Döngü Islak Laboratuvar (Self-Driving Lab Cycle)
ax_left.set_xlim(-2.5, 2.5)
ax_left.set_ylim(-2.5, 2.5)
ax_left.axis('off')
ax_left.set_title('Otonom Kapalı Döngü Islak Laboratuvar (Self-Driving Lab)', color='#ffffff', fontsize=12, fontweight='bold', pad=15)

# Döngü merkez çemberi
circle_bg = Circle((0, 0), 1.6, fill=False, edgecolor='#334155', linestyle='--', linewidth=1.5)
ax_left.add_patch(circle_bg)

loop_nodes = [
    (0, 1.6, 'AI Agent Hipotezi\n(LLM & Difuzyon)', '#38bdf8'),
    (1.6, 0, 'Robotik Sivi Isleme\n& 4D Biyobaski', '#10b981'),
    (0, -1.6, 'Otomatik scRNA-seq\n& Floresan Tarama', '#f59e0b'),
    (-1.6, 0, 'Aktif Ogrenme ile\nProtokol Guncellemesi', '#ec4899'),
]

for nx, ny, text, col in loop_nodes:
    box = FancyBboxPatch((nx - 0.75, ny - 0.35), 1.5, 0.7, boxstyle="round,pad=0.08", fc='#1e293b', ec=col, lw=1.8)
    ax_left.add_patch(box)
    ax_left.text(nx, ny, text, color='#f8fafc', fontsize=7.5, fontweight='bold', ha='center', va='center')

# Oklar
angles = [np.pi/2, 0, -np.pi/2, np.pi]
for ang in angles:
    rad1 = ang - 0.35
    rad2 = ang - 1.2
    x1, y1 = 1.6 * np.cos(rad1), 1.6 * np.sin(rad1)
    x2, y2 = 1.6 * np.cos(rad2), 1.6 * np.sin(rad2)
    arrow = FancyArrowPatch((x1, y1), (x2, y2), connectionstyle=f"arc3,rad=0.2", arrowstyle='->', mutation_scale=16, color='#60a5fa', lw=2)
    ax_left.add_patch(arrow)

ax_left.text(0, 0, 'Gunde 500+\nOtonom Iterasyon', color='#94a3b8', fontsize=9, fontweight='bold', ha='center', va='center')

# Sağ Taraf: In Vivo Epigenetik Rejuvenasyon (Horvath Saati Geri Sayımı)
chrono_age = np.linspace(20, 80, 60)
# Normal yaşlanma
bio_age_normal = chrono_age + np.random.normal(0, 1.2, len(chrono_age))
# AI-Güdümlü Periyodik Pulsed OSKM Müdahalesi
bio_age_rejuv = np.copy(bio_age_normal)
mask_intervene = chrono_age >= 50
bio_age_rejuv[mask_intervene] = 50 - 0.45 * (chrono_age[mask_intervene] - 50) + np.random.normal(0, 1.0, np.sum(mask_intervene))

ax_right.plot(chrono_age, bio_age_normal, color='#ef4444', linestyle='--', linewidth=2.0, label='Müdahalesiz Normal Epigenetik Yaşlanma')
ax_right.plot(chrono_age, bio_age_rejuv, color='#10b981', linewidth=3.0, label='In Vivo Pulsed OSKM ile Epigenetik Gençleşme')
ax_right.axvline(50, color='#f59e0b', linestyle=':', linewidth=1.8, label='Periyodik AI Rejuvenasyon Başlangıcı')

ax_right.set_title('In Vivo Epigenetik Saatin (Horvath Clock) Geri Çevrilmesi', color='#ffffff', fontsize=12, fontweight='bold', pad=15)
ax_right.set_xlabel('Kronolojik Yaş (Yıl)', color='#94a3b8', fontsize=10)
ax_right.set_ylabel('Biyolojik / Epigenetik Yaş (Yıl)', color='#94a3b8', fontsize=10)
ax_right.tick_params(colors='#94a3b8')
for spine in ax_right.spines.values():
    spine.set_color('#334155')

ax_right.grid(True, linestyle='--', alpha=0.25, color='#475569')
ax_right.legend(loc='upper left', frameon=True, facecolor='#1e293b', edgecolor='#334155', labelcolor='#e2e8f0', fontsize=8.5)

# Açıklama kutucuğu
rejuv_box = FancyBboxPatch((52, 24), 26, 12, boxstyle="round,pad=0.5", fc='#1e293b', ec='#10b981', lw=1.2)
ax_right.add_patch(rejuv_box)
ax_right.text(65, 30, 'Hucre Kimligi Korunur\n(Teratom Yok, Sadece Genclesme)\nKronolojik 75 Yas -> Biyolojik 38 Yas', color='#f8fafc', fontsize=8, fontweight='bold', ha='center', va='center')

plt.suptitle('Rejeneratif Tıbbın Nihai Sınırı: Otonom Keşif ve Hücresel Gençleşme', color='#ffffff', fontsize=14, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig('assets/06_closed_loop_self_driving_lab.png', dpi=300, facecolor=fig.get_facecolor(), edgecolor='none', bbox_inches='tight')
plt.close()
print('Vision 2 created: assets/06_closed_loop_self_driving_lab.png')
