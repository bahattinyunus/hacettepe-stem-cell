"""
🧬 Hacettepe Kök Hücre ve Rejeneratif Tıp İnteraktif Araştırma Platformu
Hacettepe Üniversitesi Kök Hücre Bilimleri Çalıştayı Notları ve Genişletilmiş Biyoinformatik Deposu.
"""

import sys
import os
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

# Proje ana dizinini sys.path'e ekle
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from src.waddington import WaddingtonLandscape, simulate_reprogramming
from src.scrna_pipeline import (
    ScRNAPipeline,
    generate_synthetic_stem_cell_data,
    CURATED_MARKER_GENES,
)
from src.stemness_index import calculate_stemness_index
from src.biomarker_discovery import DifferentialExpression, BiomarkerClassifier
from src.secretome_network import SecretomeNetwork, PARACRINE_DATABASE

# Sayfa yapılandırması
st.set_page_config(
    page_title="Hacettepe Kök Hücre & Rejeneratif Tıp Platformu",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Özel CSS stilleri
st.markdown(
    """
    <style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1E3A8A;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 1.1rem;
        color: #475569;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background-color: #F8FAFC;
        border-radius: 8px;
        padding: 1rem;
        border: 1px solid #E2E8F0;
    }
    .citation-box {
        background-color: #EFF6FF;
        border-left: 4px solid #3B82F6;
        padding: 0.8rem 1rem;
        border-radius: 4px;
        margin-bottom: 1rem;
        font-size: 0.95rem;
    }
    </style>
    """,
    unsafe_allow_html=True,
)


@st.cache_data
def get_cached_scrna_data():
    """Örnek tek hücre veri kümesini önbelleğe alır."""
    counts_df, metadata_df = generate_synthetic_stem_cell_data(n_cells=650, seed=42)
    pipeline = ScRNAPipeline()
    log1p_df = pipeline.preprocess(counts_df)
    pipeline.find_highly_variable_genes(top_n=30)
    pca_df = pipeline.run_pca(n_components=5)
    clusters = pipeline.cluster_cells(n_clusters=6)
    pseudotime = pipeline.compute_pseudotime()
    stemness = calculate_stemness_index(log1p_df)

    metadata_df["PC_1"] = pca_df["PC_1"]
    metadata_df["PC_2"] = pca_df["PC_2"]
    metadata_df["PC_3"] = pca_df["PC_3"]
    metadata_df["cluster"] = clusters.astype(str)
    metadata_df["pseudotime"] = pseudotime
    metadata_df["stemness_index"] = stemness

    return counts_df, metadata_df, log1p_df, pipeline


def main():
    st.sidebar.image("https://img.icons8.com/color/96/dna-helix.png", width=70)
    st.sidebar.title("🧬 Hacettepe Kök Hücre")
    st.sidebar.caption("Rejeneratif Tıp & Biyoinformatik Çalıştayı")

    menu = st.sidebar.radio(
        "Modül ve Analiz Seçimi:",
        [
            "🏠 Genel Bakış & Ekosistem",
            "🏔️ Waddington & Epigenetik Reprogramlama",
            "🔬 scRNA-seq & Hücre Kaderi Analizi",
            "🦠 Kanser Kök Hücreleri & THY1/CD90",
            "🧪 Rejeneratif Kaos & MSC Sekretom Ağı",
            "📊 Multi-Omik & Biyobelirteç Keşfi",
            "🎓 Çalıştay Sınavı & Bilgi Kartları",
        ],
    )

    st.sidebar.divider()
    st.sidebar.info(
        "**Hacettepe Üniversitesi**\n"
        "Kök Hücre Bilimleri Ana Bilim Dalı\n"
        "Çalıştay Araştırma ve Biyoinformatik Modülü"
    )

    counts_df, metadata_df, log1p_df, pipeline = get_cached_scrna_data()

    # ---------------------------------------------------------
    # 1. GENEL BAKIŞ & EKOSİSTEM
    # ---------------------------------------------------------
    if menu == "🏠 Genel Bakış & Ekosistem":
        st.markdown('<div class="main-header">🧬 Hacettepe Kök Hücre ve Rejeneratif Tıp Bilgi Bankası</div>', unsafe_allow_html=True)
        st.markdown('<div class="sub-header">Biyoloji, Epigenetik, Nanotıp ve Yapay Zekayı Buluşturan Disiplinler Arası Platform</div>', unsafe_allow_html=True)

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Toplam Analiz Modülü", "5 Modül", "Temelden İleriye")
        col2.metric("Tek Hücre (scRNA)", f"{len(metadata_df)} Hücre", "7 Fenotip")
        col3.metric("Biyobelirteç Paneli", f"{counts_df.shape[1]} Gen", "Multi-Omik")
        col4.metric("Parakrin Etkileşim", f"{len(PARACRINE_DATABASE)} Aks", "MSC Sekretom")

        st.divider()

        st.markdown("### 📋 Platformun Kapsamı ve Hacettepe Çalıştayı Oturumları")

        col_a, col_b = st.columns(2)
        with col_a:
            st.markdown(
                """
                #### 🧬 I. Oturum: Hücre Kimliğini Yeniden Yazmak ve İmmün Stratejiler
                * **Yeniden Programlamada Epigenetiğin Rolü** (*Öğr. Gör. Dr. Burcu Özçimen Batman*): Somatik hücrelerin epigenetik bariyeri aşarak iPSC'ye dönüşümü.
                * **İmmün Yetmezlikler için Gen Tedavisi** (*Dr. Özgür Doğuş Erol Duyar*): Ex-vivo gen tedavisi ve kliniğe translasyon.
                * **Hücre Göçü ve Otoinflamasyon** (*Dr. Öğr. Üyesi Tayfun Hilmi Akbaba*): Hastalıkların in vitro modellenmesi.

                #### 🔬 II. Oturum: Rejeneratif Kaos ve Mikroçevre
                * **Podosit Hasarında MSC Rejeneratif Etkisi** (*Öğr. Gör. Dr. Özge Burcu Şahan*): Böbrek glomerülünde parakrin ve eksozom kurtarması.
                * **Kemik İliği Stromasında Genotoksik Stres** (*Dr. Bihter Muratoğlu*): Hücre kaderi ve mikroçevre hasarı.
                * **Deneysel Pulmoner Fibroziste Yeni Yaklaşımlar** (*Dr. Öğr. Üyesi Yasemin Karaman*).
                """
            )

        with col_b:
            st.markdown(
                """
                #### 🦠 III. Oturum: Biyobelirteçler ve Terapötik Hedefler
                * **Kanser Kök Hücrelerinin Hedeflenmesi: Biyomimetik Yaklaşımlar** (*Araş. Gör. Dr. Özlem Altundağ Erdoğan*): Akıllı nanopartiküllerle CSC yok edilmesi.
                * **Kolorektal Kanserde THY1/CD90** (*Araş. Gör. Dr. Kemal Kayaşlı*): Tümör mikroçevresi, kemorezistans ve CSC kimliği.
                * **Multi-Omik Entegrasyon: Binlerce Sinyal, Kaç Biyobelirteç?** (*Dr. Öğr. Üyesi Cemil Can Eylem*): Büyük veriden hassas tıbba.

                #### 🌐 IV. Oturum: AB Fonları & Gelecek Vizyonu
                * **ZEBRARE-TWIN** (*Prof. Dr. Pervin Rukiye Dinçer*): Nadir hastalıklarda zebra balığı modelleri.
                * **OPENRAM** (*Prof. Dr. Zehra Taşkın*) & **BEEXPERT** (*Araş. Gör. Dr. Aygün Schiesser*).
                * **AI Çağında Akademi**: Islak laboratuvar (wet-lab) ile kuru laboratuvar (dry-lab) sinerjisi.
                """
            )

        st.info("Sol menüden dilediğiniz analiz modülünü seçerek etkileşimli deneyleri başlatabilirsiniz.")

        st.divider()
        st.markdown("### 🔮 AI Çağında Kök Hücre ve Rejeneratif Tıp: 2030 ve Ötesi Vizyonu")
        v_col1, v_col2, v_col3 = st.columns(3)
        with v_col1:
            st.markdown(
                """
                **1. Hücrenin 'Dijital İkizi'**  
                scGPT, Geneformer ve AlphaFold 3 entegrasyonu ile hastaya özel hücrelerin in silico modellenmesi ve reprogramlama yörüngelerinin saniyeler içinde simülasyonu.
                
                **2. De Novo Organogenez & 4D Biyobaskı**  
                Üretken yapay zeka ile kılcal damar ağları optimize edilmiş canlı biyoböbrek ve karaciğer dokusu üretimi.
                """
            )
        with v_col2:
            st.markdown(
                """
                **3. Otonom Islak Laboratuvarlar (Self-Driving Labs)**  
                AI ajanlarının hipotez kurup robotik pipetleme ve scRNA-seq döngüleriyle kapalı döngüde (closed-loop) protokolleri optimize etmesi.

                **4. Kanser Kök Hücrelerinde Diferansiasyon Terapisi**  
                CSC'leri toksik kemoterapiyle bombalamak yerine, Waddington yüzeyinde epigenetik anahtarlarla zararsız somatik hücreye dönüştürme.
                """
            )
        with v_col3:
            st.markdown(
                """
                **5. In Vivo Epigenetik Rejuvenasyon**  
                Yamanaka faktörlerinin (OSKM) geçici pulsatil ifadesiyle hücresel kimlik korunarak Horvath epigenetik saatinin 20-30 yıl gençleştirilmesi.

                **6. Sentetik Embriyoloji & AI-Biyoetik**  
                Blastoid ve gastruloid araştırmalarında uluslararası açık bilim (OPENRAM) standartları ve küresel sağlık eşitliği.
                """
            )

    # ---------------------------------------------------------
    # 2. WADDINGTON & EPİGENETİK REPROGRAMLAMA
    # ---------------------------------------------------------
    elif menu == "🏔️ Waddington & Epigenetik Reprogramlama":
        st.header("🏔️ Conrad Waddington Epigenetik Manzarası ve Reprogramlama")
        st.markdown(
            '<div class="citation-box"><b>Oturum I Araştırma Notu:</b> Öğr. Gör. Dr. Burcu Özçimen Batman sunumunda '
            'hücre kimliğinin geri döndürülebilir olduğunu ve Yamanaka Faktörleri (Oct4, Sox2, Klf4, c-Myc) ile epigenetik '
            'vadilerin aşılabildiğini vurgulamıştır.</div>',
            unsafe_allow_html=True,
        )

        col_ctrl1, col_ctrl2, col_ctrl3 = st.columns(3)
        yamanaka_input = col_ctrl1.slider("Yamanaka Kokteyli Gücü (Faktör Dozu)", 0.0, 1.0, 0.75, step=0.05)
        noise_input = col_ctrl2.slider("Hücresel Stokastisite (Termal Gürültü)", 0.02, 0.30, 0.12, step=0.02)
        steps_input = col_ctrl3.slider("Simülasyon Adım Sayısı", 50, 400, 200, step=25)

        # 3D Yüzey hesabı
        landscape = WaddingtonLandscape()
        xx, yy, zz = landscape.generate_surface_grid(
            grid_points=40,
            yamanaka_factor_strength=yamanaka_input,
        )

        # Reprogramlama simülasyonu
        sim_res = simulate_reprogramming(
            start_pos=(1.4, 0.7),
            yamanaka_strength=yamanaka_input,
            steps=steps_input,
            diffusion_noise=noise_input,
            seed=42,
        )
        traj = sim_res["trajectory"]
        traj_z = landscape.potential(traj[:, 0], traj[:, 1], yamanaka_factor_strength=yamanaka_input)

        # 3D Plotly Grafiği
        fig_3d = go.Figure()
        fig_3d.add_trace(
            go.Surface(
                x=xx,
                y=yy,
                z=zz,
                colorscale="Viridis",
                opacity=0.85,
                name="Epigenetik Potansiyel Enerji",
                colorbar=dict(title="V(x, y)"),
            )
        )
        fig_3d.add_trace(
            go.Scatter3d(
                x=traj[:, 0],
                y=traj[:, 1],
                z=traj_z + 0.15,
                mode="lines+markers",
                marker=dict(size=3, color="red"),
                line=dict(color="red", width=5),
                name="Hücre Yörüngesi (Cell Trajectory)",
            )
        )
        # Başlangıç ve bitiş noktaları
        fig_3d.add_trace(
            go.Scatter3d(
                x=[traj[0, 0]],
                y=[traj[0, 1]],
                z=[traj_z[0] + 0.3],
                mode="markers+text",
                marker=dict(size=8, color="blue"),
                text=["Somatik Hücre"],
                textposition="top center",
                name="Başlangıç (Fibroblast)",
            )
        )
        fig_3d.add_trace(
            go.Scatter3d(
                x=[traj[-1, 0]],
                y=[traj[-1, 1]],
                z=[traj_z[-1] + 0.3],
                mode="markers+text",
                marker=dict(size=8, color="green" if sim_res["success"] else "orange"),
                text=[sim_res["final_state"]],
                textposition="top center",
                name="Son Durum",
            )
        )

        fig_3d.update_layout(
            title="Waddington Potansiyel Enerji Manzarası ve Reprogramlama Dinamiği",
            scene=dict(
                xaxis_title="Farklılaşma Ekseni (<- Pluripotensi | Somatik ->)",
                yaxis_title="Doku / Soy Dallanması",
                zaxis_title="Epigenetik Potansiyel V(x, y)",
            ),
            height=600,
            margin=dict(l=0, r=0, b=0, t=40),
        )

        st.plotly_chart(fig_3d, use_container_width=True)

        res_col1, res_col2 = st.columns(2)
        with res_col1:
            st.metric("Son Hücre Durumu", sim_res["final_state"])
            if sim_res["success"]:
                st.success("✅ Tebrikler! Yamanaka faktörleri epigenetik bariyeri düşürdü ve hücre pluripotensiye kavuştu (iPSC oluştu).")
            else:
                st.warning("⚠️ Yamanaka dozu yetersiz kaldı veya hücre epigenetik vadiyi aşamayarak somatik/kısmi durumda takılı kaldı.")

        with res_col2:
            fig_energy = px.line(
                x=sim_res["steps"],
                y=sim_res["energy"],
                labels={"x": "Zaman Adımı", "y": "Potansiyel Enerji"},
                title="Hücrenin Epigenetik Enerji Seyri",
            )
            st.plotly_chart(fig_energy, use_container_width=True)

    # ---------------------------------------------------------
    # 3. scRNA-seq & HÜCRE KADERİ ANALİZİ
    # ---------------------------------------------------------
    elif menu == "🔬 scRNA-seq & Hücre Kaderi Analizi":
        st.header("🔬 Tek Hücre RNA Dizileme (scRNA-seq) Laboratuvarı")
        st.markdown(
            '<div class="citation-box">Farklılaşma ve reprogramlama süreçlerinin moleküler çözünürlüğü tek hücre '
            'transkriptomik ile aydınlatılır. Hücreler UMAP/PCA uzayında kümelenir ve gen ekspresyon düzeyleri incelenir.</div>',
            unsafe_allow_html=True,
        )

        tab1, tab2, tab3 = st.tabs(["📊 2D/3D Hücre Haritası", "🧬 Gen Ekspresyon Dağılımı", "⏱️ Pseudotime Trajectory"])

        with tab1:
            col_dim1, col_dim2 = st.columns([1, 3])
            with col_dim1:
                color_by = st.selectbox(
                    "Hücreleri Renklendir:",
                    ["cell_type", "cluster", "stemness_index", "pseudotime", "mito_percent"],
                )
                dim_type = st.radio("Boyut:", ["2D (PC1 - PC2)", "3D (PC1 - PC2 - PC3)"])

            with col_dim2:
                if "2D" in dim_type:
                    fig_pca = px.scatter(
                        metadata_df,
                        x="PC_1",
                        y="PC_2",
                        color=color_by,
                        hover_data=["cell_id", "cell_type", "total_counts", "detected_genes"],
                        title=f"scRNA-seq İzdüşümü ({color_by} ile renklendirilmiş)",
                        color_continuous_scale="Viridis",
                        height=550,
                    )
                else:
                    fig_pca = px.scatter_3d(
                        metadata_df,
                        x="PC_1",
                        y="PC_2",
                        z="PC_3",
                        color=color_by,
                        hover_data=["cell_id", "cell_type"],
                        title=f"3D scRNA-seq İzdüşümü ({color_by})",
                        color_continuous_scale="Viridis",
                        height=600,
                    )
                st.plotly_chart(fig_pca, use_container_width=True)

        with tab2:
            st.subheader("Hücre Tiplerine Göre Belirteç Gen Ekspresyonları")
            all_genes = list(log1p_df.columns)
            selected_gene = st.selectbox("İncelenecek Geni Seçin:", all_genes, index=all_genes.index("POU5F1") if "POU5F1" in all_genes else 0)

            plot_df = metadata_df.copy()
            plot_df["expression"] = log1p_df[selected_gene]

            col_v1, col_v2 = st.columns(2)
            with col_v1:
                fig_box = px.violin(
                    plot_df,
                    x="cell_type",
                    y="expression",
                    color="cell_type",
                    box=True,
                    points="all",
                    title=f"{selected_gene} Geninin Fenotipler Arası Ekspresyon Dağılımı",
                    labels={"expression": "Log1p Normalize Ekspresyon", "cell_type": "Hücre Tipi"},
                )
                fig_box.update_layout(showlegend=False, xaxis_tickangle=-40)
                st.plotly_chart(fig_box, use_container_width=True)

            with col_v2:
                fig_feature = px.scatter(
                    plot_df,
                    x="PC_1",
                    y="PC_2",
                    color="expression",
                    color_continuous_scale="Reds",
                    title=f"Hücre Haritası Üzerinde {selected_gene} Ekspresyon Yoğunluğu",
                )
                st.plotly_chart(fig_feature, use_container_width=True)

        with tab3:
            st.subheader("⏱️ Pseudotime: Kök Hücreden Farklılaşmaya Gelişim Zamanı")
            fig_pt = px.scatter(
                metadata_df,
                x="pseudotime",
                y="stemness_index",
                color="cell_type",
                hover_data=["cell_id"],
                title="Hücresel Gelişim Zamanı (Pseudotime) ile Kök Hücrelilik İndeksi (mRNAsI) İlişkisi",
                labels={"pseudotime": "Pseudotime (0: Kök Hücre Kökü -> 1: Nihai Farklılaşmış)", "stemness_index": "mRNAsI Kök Hücrelilik"},
            )
            st.plotly_chart(fig_pt, use_container_width=True)
            st.caption("Grafik, hücrelerin diferansiasyon geçirdikçe pluripotensi skorunun azaldığını doğrular.")

    # ---------------------------------------------------------
    # 4. KANSER KÖK HÜCRELERİ & THY1/CD90
    # ---------------------------------------------------------
    elif menu == "🦠 Kanser Kök Hücreleri & THY1/CD90":
        st.header("🦠 Kanser Kök Hücreleri (CSC) ve THY1/CD90 Hedefleme")
        st.markdown(
            '<div class="citation-box"><b>Oturum III Notu:</b> Araş. Gör. Dr. Kemal Kayaşlı, kolorektal kanserde THY1/CD90 '
            'ekspresyonunun kanser kök hücresi fenotipi, tümör nüksü ve kemorezistans ile doğrudan bağlantılı olduğunu göstermiştir. '
            'Araş. Gör. Dr. Özlem Altundağ Erdoğan ise biyomimetik taşıyıcılarla bu hücrelerin hedeflenmesini sunmuştur.</div>',
            unsafe_allow_html=True,
        )

        col_c1, col_c2 = st.columns([1, 2])
        with col_c1:
            st.subheader("⚙️ Hedefli Tedavi Parametreleri")
            thy1_threshold = st.slider("THY1/CD90 Pozitiflik Eşik Değeri", 1.0, 5.0, 2.8, step=0.2)
            nanoparticle_targeting = st.slider("Biyomimetik Nanopartikül Hedefleme Gücü (%)", 50, 99, 88)
            chemo_drug = st.selectbox("Kemoterapi Ajanı:", ["5-Fluorouracil (5-FU)", "Oxaliplatin", "Irinotecan"])

        # CSC Analiz Verisi
        crc_mask = metadata_df["cell_type"].isin(["Kolorektal CSC", "Somatik Fibroblast", "Hasarli Podosit"])
        crc_cells = metadata_df[crc_mask].copy()
        crc_cells["THY1_expr"] = log1p_df.loc[crc_cells.index, "THY1"]
        crc_cells["CD90_Positive"] = crc_cells["THY1_expr"] >= thy1_threshold

        # Kemorezistans skoru simülasyonu
        crc_cells["chemo_resistance"] = (crc_cells["THY1_expr"] * 0.45 + crc_cells["stemness_index"] * 0.55).clip(0, 1)

        with col_c2:
            fig_resist = px.scatter(
                crc_cells,
                x="THY1_expr",
                y="chemo_resistance",
                color="cell_type",
                size="stemness_index",
                labels={"THY1_expr": "THY1/CD90 Log1p Ekspresyonu", "chemo_resistance": f"{chemo_drug} Direnç Skoru"},
                title=f"THY1/CD90 Ekspresyonu ile {chemo_drug} Direnci Arasındaki Korelasyon",
            )
            fig_resist.add_vline(x=thy1_threshold, line_dash="dash", line_color="red", annotation_text="CSC Pozitif Eşik")
            st.plotly_chart(fig_resist, use_container_width=True)

        st.subheader("🎯 Biyomimetik Taşıyıcı Terapötik Yanıt Simülasyonu")
        col_res1, col_res2, col_res3 = st.columns(3)
        csc_count = int(crc_cells["CD90_Positive"].sum())
        eradicated_csc = int(csc_count * (nanoparticle_targeting / 100.0))
        recurrence_risk = max(5, int((1.0 - (eradicated_csc / max(1, csc_count))) * 100))

        col_res1.metric("Toplam Tespit Edilen CSC (THY1+)", f"{csc_count} Hücre")
        col_res2.metric("Yok Edilen Kanser Kök Hücresi", f"{eradicated_csc} Hücre", f"%{nanoparticle_targeting} etkinlik")
        col_res3.metric("Tahmini Tümör Nüks Riski", f"%{recurrence_risk}", delta=f"-{nanoparticle_targeting}%", delta_color="inverse")

    # ---------------------------------------------------------
    # 5. REJENERATİF KAOS & MSC SEKRETOM AĞI
    # ---------------------------------------------------------
    elif menu == "🧪 Rejeneratif Kaos & MSC Sekretom Ağı":
        st.header("🧪 Rejeneratif Kaos, Podosit Hasarı & MSC Sekretom Ağı")
        st.markdown(
            '<div class="citation-box"><b>Oturum II Notu:</b> Öğr. Gör. Dr. Özge Burcu Şahan sunumunda '
            'kemik iliği kaynaklı mezenkimal kök hücrelerin (BM-MSC) ve salgıladıkları eksozomların '
            'podosit hasarı mikroçevresindeki onarıcı rolünü açıklamıştır.</div>',
            unsafe_allow_html=True,
        )

        sec_net = SecretomeNetwork()
        secretome_df = sec_net.get_dataframe()

        col_p1, col_p2 = st.columns([1, 2])
        with col_p1:
            st.subheader("🛠️ Tedavi Parametreleri")
            msc_dose = st.slider("MSC Hücre Dozu Katsayısı", 0.5, 2.5, 1.2, step=0.1)
            exosome_enriched = st.checkbox("Saflaştırılmış Eksozom Takviyesi Ekle", value=True)
            initial_damage = st.slider("Glomerüler Hasar / Proteinüri Düzeyi", 0.3, 1.0, 0.85, step=0.05)

            rescue_sim = sec_net.simulate_podocyte_rescue(
                msc_dose=msc_dose,
                exosome_enriched=exosome_enriched,
                initial_damage_score=initial_damage,
            )

        with col_p2:
            timeline = rescue_sim["timeline"]
            fig_time = go.Figure()
            fig_time.add_trace(
                go.Scatter(
                    x=timeline["days"],
                    y=timeline["damage_curve"],
                    mode="lines",
                    name="Doku Hasarı / Kaos Skoru",
                    line=dict(color="red", width=3),
                )
            )
            fig_time.add_trace(
                go.Scatter(
                    x=timeline["days"],
                    y=timeline["nephrin_curve"],
                    mode="lines",
                    name="Nefrin (NPHS1) İyileşme Skoru",
                    line=dict(color="green", width=3),
                )
            )
            fig_time.update_layout(
                title="28 Günlük Tedavi Süresince Podosit İyileşme Eğrisi",
                xaxis_title="Gün",
                yaxis_title="Skor (0 - 1)",
            )
            st.plotly_chart(fig_time, use_container_width=True)

        m1, m2, m3 = st.columns(3)
        m1.metric("Proteinüri Düşüş Oranı", f"%{rescue_sim['proteinuria_reduction_pct']:.1f}")
        m2.metric("Nefrin (NPHS1) Restorasyonu", f"%{rescue_sim['nephrin_recovery']*100:.1f}")
        m3.metric("Kalan Mikroçevre Hasarı", f"{rescue_sim['final_damage']:.2f}")

        st.subheader("🔗 MSC Parakrin ve Eksozom Etkileşim Haritası")
        st.dataframe(
            secretome_df[["ligand", "receptor", "target_cell", "pathway", "action", "organ"]],
            use_container_width=True,
        )

    # ---------------------------------------------------------
    # 6. MULTİ-OMİK & BİYOBELİRTEÇ KEŞFİ
    # ---------------------------------------------------------
    elif menu == "📊 Multi-Omik & Biyobelirteç Keşfi":
        st.header("📊 Multi-Omik Entegrasyon ve Biyobelirteç Keşfi (Dry-Lab)")
        st.markdown(
            '<div class="citation-box"><b>Oturum III Araştırma Notu:</b> Dr. Öğr. Üyesi Cemil Can Eylem, '
            '"Binlerce Sinyal, Kaç Biyobelirteç?" başlığıyla omik verilerinden (transkriptom, proteom) klinik olarak '
            'anlamlı minimal biyobelirteç setlerinin yapay zeka ile nasıl süzülebileceğini açıklamıştır.</div>',
            unsafe_allow_html=True,
        )

        tab_de, tab_ml = st.tabs(["🌋 Diferansiyel Ekspresyon & Volcano Plot", "🤖 Makine Öğrenmesi ile Biyobelirteç Seçimi"])

        with tab_de:
            col_g1, col_g2 = st.columns(2)
            all_types = metadata_df["cell_type"].unique().tolist()
            grp_a_type = col_g1.selectbox("Deney Grubu:", all_types, index=all_types.index("Kolorektal CSC") if "Kolorektal CSC" in all_types else 0)
            grp_b_type = col_g2.selectbox("Kontrol Grubu:", all_types, index=all_types.index("Somatik Fibroblast") if "Somatik Fibroblast" in all_types else 1)

            grp_a_cells = metadata_df[metadata_df["cell_type"] == grp_a_type].index.tolist()
            grp_b_cells = metadata_df[metadata_df["cell_type"] == grp_b_type].index.tolist()

            col_th1, col_th2 = st.columns(2)
            fc_th = col_th1.slider("Log2 Fold Change Eşiği", 0.5, 3.0, 1.2, step=0.1)
            q_th = col_th2.slider("FDR Q-değeri Eşiği", 0.001, 0.1, 0.05, step=0.005)

            de_df = DifferentialExpression.analyze(
                log1p_df,
                group_a_samples=grp_a_cells,
                group_b_samples=grp_b_cells,
                log2_fc_threshold=fc_th,
                fdr_threshold=q_th,
            )

            fig_volcano = px.scatter(
                de_df,
                x="log2_fold_change",
                y="-log10_p_value",
                color="status",
                hover_data=["gene", "fdr_q_value"],
                color_discrete_map={
                    "Up-regulated": "#EF4444",
                    "Down-regulated": "#3B82F6",
                    "Not Significant": "#94A3B8",
                },
                title=f"Volcano Plot: {grp_a_type} vs {grp_b_type}",
                labels={"log2_fold_change": "Log2 Fold Change", "-log10_p_value": "-Log10 P-değeri"},
                height=520,
            )
            fig_volcano.add_vline(x=fc_th, line_dash="dash", line_color="gray")
            fig_volcano.add_vline(x=-fc_th, line_dash="dash", line_color="gray")
            st.plotly_chart(fig_volcano, use_container_width=True)

            st.dataframe(
                de_df[de_df["status"] != "Not Significant"][["gene", "log2_fold_change", "fdr_q_value", "status"]],
                use_container_width=True,
            )

        with tab_ml:
            st.subheader("🤖 Random Forest & Özellik Önemi ile Minimal Biyobelirteç Paneli")
            
            # CSC vs Somatik sınıflandırma görevi
            task_mask = metadata_df["cell_type"].isin([grp_a_type, grp_b_type])
            sub_meta = metadata_df[task_mask]
            X_task = log1p_df.loc[sub_meta.index]
            y_task = (sub_meta["cell_type"] == grp_a_type).astype(int)

            top_k_input = st.slider("Seçilecek Biyobelirteç Sayısı (Minimal Panel)", 3, 15, 6)
            clf = BiomarkerClassifier(n_estimators=100, random_state=42)
            importances = clf.fit_and_rank(X_task, y_task, top_k=top_k_input)
            cv_eval = clf.evaluate_cv(X_task, y_task, n_splits=5)

            col_m1, col_m2 = st.columns(2)
            with col_m1:
                st.metric("Çapraz Doğrulama ROC-AUC", f"{cv_eval['roc_auc']:.3f}")
                st.metric("Model Doğruluğu (Accuracy)", f"%{cv_eval['accuracy']*100:.1f}")

                fig_roc = px.area(
                    x=cv_eval["fpr"],
                    y=cv_eval["tpr"],
                    title=f"Biyobelirteç Paneli ROC Eğrisi (AUC = {cv_eval['roc_auc']:.3f})",
                    labels={"x": "Yanlış Pozitif Oranı (1 - Özgüllük)", "y": "Doğru Pozitif Oranı (Duyarlılık)"},
                )
                fig_roc.add_shape(type="line", line=dict(dash="dash"), x0=0, x1=1, y0=0, y1=1)
                st.plotly_chart(fig_roc, use_container_width=True)

            with col_m2:
                top_genes_df = pd.DataFrame({
                    "Gen": importances.head(top_k_input).index,
                    "Gini Önemi": importances.head(top_k_input).values,
                })
                fig_imp = px.bar(
                    top_genes_df,
                    x="Gini Önemi",
                    y="Gen",
                    orientation="h",
                    title=f"En Ayırt Edici {top_k_input} Biyobelirteç",
                    color="Gini Önemi",
                    color_continuous_scale="Blues",
                )
                fig_imp.update_layout(yaxis=dict(autorange="reversed"))
                st.plotly_chart(fig_imp, use_container_width=True)

    # ---------------------------------------------------------
    # 7. ÇALIŞTAY SINAVI & BİLGİ KARTLARI
    # ---------------------------------------------------------
    elif menu == "🎓 Çalıştay Sınavı & Bilgi Kartları":
        st.header("🎓 Hacettepe Kök Hücre Çalıştayı İnteraktif Sınavı")
        st.markdown("Bilgilerinizi çalıştayın 4 oturumundaki temel kazanımlara göre test edin.")

        questions = [
            {
                "q": "1. 2012 Nobel Tıp Ödülü'nü kazandıran ve somatik hücreleri pluripotent duruma (iPSC) geriye programlayan Yamanaka Faktörleri hangileridir?",
                "options": [
                    "Oct3/4, Sox2, Klf4, c-Myc",
                    "p53, Rb, BRCA1, APC",
                    "CD73, CD90, CD105, CD44",
                    "VEGF, HGF, TGFb, IL10",
                ],
                "answer": "Oct3/4, Sox2, Klf4, c-Myc",
                "explain": "Shinya Yamanaka, bu 4 transkripsiyon faktörünü retroviral vektörlerle somatik hücreye aktararak iPSC hücrelerini üretmiştir.",
            },
            {
                "q": "2. Conrad Waddington'ın Epigenetik Manzara modelinde, farklılaşmış bir somatik hücrenin durumu nasıl tanımlanır?",
                "options": [
                    "Vadinin en derin tabanında yer alan kararlı minimum enerji durumu",
                    "Tepenin en zirvesinde yer alan totipotent potansiyel",
                    "Sürekli hareket eden belirsiz faz",
                    "Enerji seviyesi sonsuza giden durum",
                ],
                "answer": "Vadinin en derin tabanında yer alan kararlı minimum enerji durumu",
                "explain": "Farklılaşma, hücrenin vadiden aşağı yuvarlanarak minimum enerjili kararlı somatik kuyuya yerleşmesidir.",
            },
            {
                "q": "3. Kemik iliği kaynaklı Mezenkimal Kök Hücrelerin (BM-MSC) hasarlı podositleri iyileştirmesindeki ana mekanizma nedir?",
                "options": [
                    "Parakrin sinyalleşme, eksozom salınımı ve anti-enflamatuvar mikroçevre modülasyonu",
                    "Doğrudan podosite transdiferansiyasyon geçirmek",
                    "Bağışıklık sistemini tamamen yok etmek",
                    "Hücre DNA dizilimini baştan aşağı değiştirmek",
                ],
                "answer": "Parakrin sinyalleşme, eksozom salınımı ve anti-enflamatuvar mikroçevre modülasyonu",
                "explain": "Arnold Caplan'ın vurguladığı üzere MSC'ler ortamda 'ilaç fabrikası' gibi davranarak salgıladıkları sekretom/eksozom ile dokuyu kurtarır.",
            },
            {
                "q": "4. Kolorektal kanserde kanser kök hücrelerini (CSC) belirlemede ve kemorezistansı anlamada önemli bir yüzey biyobelirteci nedir?",
                "options": [
                    "THY1 (CD90)",
                    "Albumin",
                    "Hemoglobin",
                    "Troponin",
                ],
                "answer": "THY1 (CD90)",
                "explain": "Araş. Gör. Dr. Kemal Kayaşlı'nın sunumunda açıklandığı gibi THY1/CD90 kolorektal CSC'lerde belirteç olarak tümör agresifliğiyle koreledir.",
            },
            {
                "q": "5. 'Multi-omik entegrasyonda binlerce sinyal arasından kaç biyobelirteç translasyona geçer?' sorusunda minimal panel seçiminde hangi yaklaşım kullanılır?",
                "options": [
                    "Makine öğrenmesi öznitelik seçimi, çapraz doğrulama ve FDR düzeltmeli istatistik",
                    "Yalnızca p-değerine bakarak ilk geni seçmek",
                    "Tüm genleri tek tek elle incelemek",
                    "Diferansiyel ekspresyonu göz ardı etmek",
                ],
                "answer": "Makine öğrenmesi öznitelik seçimi, çapraz doğrulama ve FDR düzeltmeli istatistik",
                "explain": "Dr. Öğr. Üyesi Cemil Can Eylem'in belirttiği üzere yüksek boyutlu omik verilerinden minimal biyobelirteç panellerine ulaşmak için ML ve istatistik şarttır.",
            },
        ]

        score = 0
        with st.form("quiz_form"):
            user_answers = []
            for i, q_item in enumerate(questions):
                st.markdown(f"**{q_item['q']}**")
                ans = st.radio("Seçiminiz:", q_item["options"], key=f"q_{i}")
                user_answers.append(ans)
                st.write("")

            submitted = st.form_submit_button("Sınavı Tamamla ve Puanı Gör")

        if submitted:
            for i, q_item in enumerate(questions):
                if user_answers[i] == q_item["answer"]:
                    score += 1

            pct = (score / len(questions)) * 100
            st.markdown(f"### 🎯 Başarı Puanınız: {score} / {len(questions)} (%{pct:.0f})")

            for i, q_item in enumerate(questions):
                if user_answers[i] == q_item["answer"]:
                    st.success(f"**Soru {i+1}: Doğru!** {q_item['explain']}")
                else:
                    st.error(f"**Soru {i+1}: Yanlış.** Doğru cevap: *{q_item['answer']}*. Açıklama: {q_item['explain']}")


if __name__ == "__main__":
    main()
