# 📂 Veri Kümeleri ve Metodoloji

Bu dizin, Hacettepe Üniversitesi Kök Hücre Bilimleri Çalıştayı'nda tartışılan deneysel ve biyoinformatik senaryoları simüle eden veri setlerini ve parakrin etkileşim veritabanını barındırır.

## 📄 Dosyalar

| Dosya Adı | Boyut / Satır | Açıklama |
| :--- | :--- | :--- |
| `scrna_counts.csv` | 500 hücre x 45 gen | Tek hücre RNA dizileme (scRNA-seq) UMI ham sayı matrisi |
| `scrna_metadata.csv` | 500 hücre | Hücre kimlikleri, hücre tipleri, kütüphane büyüklüğü ve mitokondriyal oranlar |
| `paracrine_interactions.csv` | 7 aks | BM-MSC sekretom ve eksozom ligand-reseptör yolakları |

## 🧬 Hücre Fenotipleri

1. **Pluripotent ESC:** Embriyonik kök hücreler (`POU5F1+`, `SOX2+`, `NANOG+`).
2. **Erken iPSC:** Yamanaka faktörleriyle indüklenmiş, henüz tam pluripotensiye ulaşmamış ara durum (`KLF4+`, `MYC+`).
3. **Geç iPSC:** Yeniden programlanması tamamlanmış indüklenmiş pluripotent kök hücreler.
4. **Somatik Fibroblast:** Doku matrisi üreten erişkin bağ dokusu hücreleri (`COL1A1+`, `ACTA2+`, `FAP+`).
5. **Kemik İliği Mezenkimal Kök Hücresi (BM-MSC):** İmmünomodülatör ve trofik hücreler (`NT5E/CD73+`, `THY1/CD90+`, `ENG/CD105+`, `PTPRC/CD45-`).
6. **Hasarlı Podosit:** Glomerüler filtrasyon bariyer hasarı taşıyan böbrek hücreleri.
7. **Kolorektal Kanser Kök Hücresi (CRC-CSC):** Kemoterapiye dirençli ve nüks potansiyeli yüksek kök hücre benzeri tümör alt popülasyonu (`THY1+`, `CD44+`, `PROM1+`).
