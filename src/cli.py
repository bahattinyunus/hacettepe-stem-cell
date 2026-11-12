"""
Hacettepe Kök Hücre ve Rejeneratif Tıp Komut Satırı Arayüzü (CLI).

Araştırmacıların veri setlerini terminal üzerinden analiz edebilmesini,
reprogramlama simülasyonları yürütebilmesini ve Streamlit platformunu
başlatabilmesini sağlayan konsol arayüzü.
"""

import argparse
import sys
import os
import json
import pandas as pd

from src.waddington import simulate_reprogramming
from src.scrna_pipeline import ScRNAPipeline, generate_synthetic_stem_cell_data
from src.stemness_index import calculate_stemness_index
from src.biomarker_discovery import DifferentialExpression, BiomarkerClassifier
from src.epigenetic_clock import HorvathEpigeneticClock
from src.pathway_enrichment import PathwayEnrichment
from src.nanomedicine_simulator import BiomimeticNanomedicine


def main():
    if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
        try:
            sys.stdout.reconfigure(encoding="utf-8")
            sys.stderr.reconfigure(encoding="utf-8")
        except Exception:
            pass

    parser = argparse.ArgumentParser(
        prog="hacettepe-stem-cell",
        description="Hacettepe Kok Hucre, Epigenetik ve Rejeneratif Tip Biyoinformatik CLI",
    )
    subparsers = parser.add_subparsers(dest="command", help="Kullanilabilir alt komutlar")

    # 1. simulate-reprogramming
    reprog_parser = subparsers.add_parser("simulate-reprogramming", help="Waddington reprogramlama simulasyonu calistirir")
    reprog_parser.add_argument("--yamanaka", type=float, default=0.85, help="Yamanaka kokteyli gucu [0-1]")
    reprog_parser.add_argument("--steps", type=int, default=150, help="Simulasyon adim sayisi")
    reprog_parser.add_argument("--noise", type=float, default=0.12, help="Termal stokastik gurultu")
    reprog_parser.add_argument("--output", type=str, default=None, help="JSON cikti dosya yolu")

    # 2. calc-stemness
    stemness_parser = subparsers.add_parser("calc-stemness", help="Ekspresyon tablosundan mRNAsI kok hucrelilik indeksi hesaplar")
    stemness_parser.add_argument("--counts", type=str, required=True, help="scRNA-seq ham sayi CSV dosyasi")
    stemness_parser.add_argument("--output", type=str, default=None, help="Sonuclarin kaydedilecegi CSV yolu")

    # 3. epigenetic-clock
    clock_parser = subparsers.add_parser("epigenetic-clock", help="Horvath epigenetik saati ve in vivo rejuvenasyon simulasyonu")
    clock_parser.add_argument("--age", type=float, default=65.0, help="Kronolojik yas")
    clock_parser.add_argument("--cycles", type=int, default=6, help="Pulsatile OSKM dongu sayisi")
    clock_parser.add_argument("--strength", type=float, default=0.85, help="Yamanaka faktor dozu")

    # 4. pathway-enrichment
    pathway_parser = subparsers.add_parser("pathway-enrichment", help="Diferansiyel genler icin biyolojik yolak zenginlesmesi (ORA) yurutur")
    pathway_parser.add_argument("--genes", type=str, required=True, help="Virgulle ayrilmis gen listesi (orn: POU5F1,SOX2,NANOG)")

    # 5. launch-app
    subparsers.add_parser("launch-app", help="Interaktif Streamlit web arayuzunu baslatir")

    args = parser.parse_args()

    if args.command == "simulate-reprogramming":
        print(f"[*] Yamanaka dozu = {args.yamanaka}, Adim = {args.steps} ile simulasyon baslatiliyor...")
        res = simulate_reprogramming(yamanaka_strength=args.yamanaka, steps=args.steps, diffusion_noise=args.noise)
        print(f"[+] Nihai Durum: {res['final_state']}")
        print(f"[+] Basari (Pluripotensiye Ulasma): {res['success']}")
        if args.output:
            data_to_save = {
                "final_state": res["final_state"],
                "success": res["success"],
                "trajectory": res["trajectory"].tolist(),
                "energy": res["energy"].tolist(),
            }
            with open(args.output, "w", encoding="utf-8") as f:
                json.dump(data_to_save, f, indent=2)
            print(f"[+] Sonuc kaydedildi: {args.output}")

    elif args.command == "calc-stemness":
        print(f"[*] Veri okunuyor: {args.counts}")
        counts = pd.read_csv(args.counts, index_col=0)
        pipeline = ScRNAPipeline()
        log1p = pipeline.preprocess(counts)
        stemness = calculate_stemness_index(log1p)
        print(f"[+] Ortalama Stemness Skoru: {stemness.mean():.4f} (Min: {stemness.min():.4f}, Max: {stemness.max():.4f})")
        if args.output:
            stemness.to_csv(args.output)
            print(f"[+] Skorlar kaydedildi: {args.output}")

    elif args.command == "epigenetic-clock":
        clock = HorvathEpigeneticClock()
        res = clock.simulate_pulsed_rejuvenation(chronological_age=args.age, cycles=args.cycles, pulse_strength=args.strength)
        print(f"[*] Kronolojik Yas: {res['chronological_age']} yil")
        print(f"[*] Baslangic Biyolojik Yasi: {res['initial_biological_age']:.1f} yil")
        print(f"[+] Rejuvenasyon Sonrasi Biyolojik Yas: {res['final_biological_age']:.1f} yil")
        print(f"[+] Kazanilan Genclesme: {res['years_rejuvenated']:.1f} yil ({res['safety_status']})")

    elif args.command == "pathway-enrichment":
        genes = [g.strip() for g in args.genes.split(",") if g.strip()]
        ora = PathwayEnrichment()
        df = ora.analyze(genes)
        if df.empty:
            print("[-] Belirtilen genler icin anlamli yolak zenginlesmesi bulunamadi.")
        else:
            print(f"[+] Zenginlesen {len(df)} Yolak:")
            print(df[["pathway", "overlap_count", "rich_factor", "p_value", "fdr_q_value"]].to_string(index=False))

    elif args.command == "launch-app":
        import subprocess
        app_path = os.path.join(os.path.dirname(__file__), "..", "app", "app.py")
        print(f"[*] Streamlit Web Platformu Baslatiliyor ({app_path})...")
        subprocess.run(["streamlit", "run", app_path])

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
