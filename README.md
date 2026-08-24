# Retrogeny człowieka - analiza genomowa

Analiza retrogenów z pliku `human_retrocopies_merged.bed` (14 874 retrogenów, hs1),
ich klasyfikacja genomowa, walidacja długimi odczytami ENCODE oraz przygotowanie
sekwencji transkryptów do analizy GC content.

## Struktura projektu

```
strict/                         ← FINALNA analiza (najważniejsza)
├── classify_4lists_strict.sh   ← klasyfikacja retrogenów na 4 listy
├── validate_encode4_strict.sh  ← walidacja ENCODE (intersekcja)
├── select_encode4_transcripts.py ← selekcja 1 transkryptu ENCODE na retrogen
├── fix_intronic_by_ensg.py     ← poprawka multi-exon (tylko intronowe)
├── list1_intergenic_strict.bed  ← 5 205 międzygenowych
├── list2_intronic_strict.bed    ← 3 191 intronowych
├── list3_cds_one_exon_strict.bed ← 768 CDS w 1 egzonie
├── list4_other_strict.bed       ← 5 710 reszta
├── list5_selfhit.bed            ← 3 256 self-hit
├── selected_encode4_{intergenic,intronic,cds}.bed ← wybrane transkrypty ENCODE
    └── fasta/
    ├── intergenic_transcripts.fa ← 941 sekwencji
    ├── intronic_transcripts.fa   ← 1 038 sekwencji
    ├── cds_transcripts.fa        ← 398 sekwencji
    ├── gc_analysis.py            ← analiza GC (okna 100 nt od 5' końca)
    ├── gc_statistics.py          ← test Mann-Whitney U
    └── gc_content_results.csv    ← surowe dane GC

ncbi_final/                     ← pośrednia wersja (elastyczne progi)
Publications/                   ← publikacje źródłowe (pliki .txt)
```

## Dane wejściowe (NIE w repo)

| Plik | Źródło |
|------|--------|
| `human_retrocopies_merged.bed` | prof. Szcześniak (CAT+Liftoff na hs1) |
| `hs1.ncbiRefSeq.bigZip.gtf.gz` | UCSC → hs1 → ncbiRefSeq |
| `hs1.fa.gz` | UCSC → hs1 → bigZips |
| `encode4_long_liftover.bed` | ENCODE4 long reads (hg38), liftover → hs1 (narzędzie `liftOver` + `hg38ToHs1.over.chain.gz`) |

## Metodyka (skrót)

1. **Klasyfikacja** - priorytet: CDS 1 egzon → intron → reszta → intergenic.
   - intronowe: 0% overlapu z eksonem (po filtrze self-overlapu po nazwie genu)
   - self-hit: retrogen mający tylko siebie jako gene_id → międzygenowe
2. **Walidacja ENCODE** - `bedtools intersect -f 0.1 -F 0.1`
3. **Selekcja transkryptów** - intronowe: single-exon > długość zbliżona > lider 5';
   intergenic/CDS: długość zbliżona > lider 5'
4. **Sekwencje** - `bedtools getfasta -split -s` (eksony poskładane, nić uwzględniona)
5. **Analiza GC** - okno 100 nt od końca 5', porównanie między kategoriami

## Wyniki analizy GC content (koniec 5')

| Kategoria | n | GC 5' (średnia) | GC całość (średnia) |
|-----------|-----|-----------------|---------------------|
| CDS       | 398 | **0.596**       | 0.524 |
| Intergenic| 941 | 0.518           | 0.476 |
| Intronic  | 1038| 0.509           | 0.468 |

**Test Mann-Whitney U (GC na 5'):**

| Porównanie | p-value | Istotność |
|------------|---------|-----------|
| Intergenic vs Intronic | 0.084 | n.s. |
| Intergenic vs CDS      | <0.001 | *** |
| Intronic vs CDS        | <0.001 | *** |

**Gradient GC wzdłuż 5'→3'** (okna 100 nt): wszystkie kategorie mają
najwyższe GC na końcu 5' i spadek w kierunku 3' — zgodne z Mordstein 2020.

## Wymagania

- bedtools, samtools
- Python 3 (skrypty selekcji)
- liftOver (tylko do regeneracji `encode4_long_liftover.bed`)
