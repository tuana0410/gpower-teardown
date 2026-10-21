# K2 Teardown Sprint - LEDGER & Log of Actions

## 1. Proje Özeti ve Kapsam
Bu proje, G*Power (v3.1.9.7) yazılımının tersine mühendislik yöntemleriyle (Ghidra analizi aracılığıyla `FUN_004e01d0` ve form başlatma rutinleri incelenerek) incelenmesi ve temel istatistiksel hesaplama motorunun sözleşmeye (`CONTRACT.md`) tam uyumlu bir CLI aracı olarak (`gpower_clone.py`) yeniden üretilmesini kapsar.

## 2. Teknik Altyapı ve Kütüphaneler
- **Dil:** Python 3.10
- **İstatistiksel Motor:** `scipy.stats` (merkezi olmayan t, F ve Chi-square dağılımları için `nct`, `ncf`, `ncx2` modülleri)
- **Arayüz / Çıktı:** Standart girdi/çıktı (CLI) üzerinden tam çift hassasiyetli JSON formatı.

## 3. Desteklenen Test Aileleri ve Modlar
- **Test Aileleri:** `t2-indep`, `t-paired`, `anova-oneway`, `chi2-contingency`, `multireg-r2`
- **Analiz Modları:** `a-priori`, `post-hoc`
- **Hata Yönetimi:** `OUT_OF_RANGE`, `MISSING_PARAM`, `INVALID_VALUE`, `NO_SOLUTION` kodları ve uygun hata çıkışları (`exit code 2` ve `1`).

## 4. Git ve Teslim Notları
Tüm simülasyon kodları, hata yakalama mekanizmaları ve parametre doğrulama kuralları `gpower_clone.py` içerisinde modüler olarak kodlanmış ve GitHub deposuna (`tuana0410/gpower-teardown`) başarıyla push edilmiştir.