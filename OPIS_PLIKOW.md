# Opis plików w archiwum

Ten plik opisuje wszystkie pliki w `analiza_retrogenow.zip`: co zawierają, do czego
służą i skąd się wzięły. Archiwum zawiera **wyniki i skrypty** analizy - nie zawiera
dużych plików wejściowych.

---

## Przebieg analizy

Pliki powstawały w tej kolejności:

```
1. dane wejściowe (retrokopie + adnotacje + odczyty ENCODE)
   ↓
2. CDS_one_exon.py          → single_cds.bed
   ↓
3. classify_4lists_strict.sh → list1–list4 + list5_selfhit.bed + unclassified_7.bed
   ↓
4. validate_encode4_strict.sh → encode4_vs_cat_list*.tsv
   ↓
5. select_encode4_transcripts.py → selected_encode4_*.bed
   ↓
6. fix_intronic_by_ensg.py      → poprawka selected_encode4_intronic.bed
   ↓
7. bedtools getfasta             → *_transcripts.fa (sekwencje)
   ↓
8. gc_analysis_halves.py         → gc_halves_results.csv
   ↓
9. gc_statistics.py              → testy statystyczne (wyniki)
   ↓
10. make_charts.py               → charts/*.png (wykresy)
```

---

## Główny poziom

### `README.md`
- opis całego projektu i metodyki.
- punkt startowy - tłumaczy strukturę, dane i kroki analizy.
- napisany ręcznie jako dokumentacja.

### `wyniki_klasyfikacji.md`
- zbiorcza tabela wyników.
- liczebność 4 kategorii (+self-hit), % wsparcia ENCODE oraz wyniki
  analizy GC (połowy, test Wilcoxona, test Mann-Whitney U z korektą Bonferroniego).
- opracowany na podstawie wyników skryptów.

### `human_retrocopies_merged.bed`
- dane wejściowe — 14 874 retrokopie.
- koordynaty genomowe (chr, start, koniec), nazwę retrogena
  (`retro_human_N|gen_rodzicielski|typ`) i nić.
- od prof. Szcześniaka (CAT + Liftoff na genomie hs1).

---

## Skrypty (`strict/`)

### `CDS_one_exon.py`
- generuje `single_cds.bed` — transkrypty z dokładnie jednym egzonem kodującym.
- czyta `ncbiRefSeq.gtf` (plik zewnętrzny, nie w ZIP).

### `classify_4lists_strict.sh`
- główna klasyfikacja retrogenów na 4 listy (+ self-hit).
- przetwarza `human_retrocopies_merged.bed` i `ncbiRefSeq.gtf`
  (bedtools intersect). Na końcu sam sprawdza, czy sumy się zgadzają.

### `validate_encode4_strict.sh`
- walidacja ekspresji — sprawdza, które retrogeny mają wsparcie
  w długich odczytach ENCODE (bedtools intersect, próg 10%).
- przecina listy retrogenów z `encode4_long_liftover.bed`.

### `select_encode4_transcripts.py`
- wybiera **jeden** transkrypt ENCODE na każdy retrogen.
- intronowe → single-exon > długość zbliżona > lider 5';
  pozostałe → długość zbliżona > lider 5'.
- czyta pliki `encode4_vs_cat_list*.tsv`.

### `fix_intronic_by_ensg.py`
- dla intronowych zamienia transkrypty multi-exon na single-exon
  z tego samego genu (ENSG), nachodzący na retrogen.
- odróżnia transkrypt retrokopii od transkryptu genu-gospodarza.
- poprawia `selected_encode4_intronic.bed`.

### `make_charts.py`
- generuje wszystkie 6 wykresów PNG.
- czyta `gc_halves_results.csv` i pliki FASTA.

---

## Wykresy (`strict/charts/`)

| Plik | Co pokazuje |
|------|-------------|
| `klasyfikacja.png` | liczba retrogenów w 3 czystych kategoriach |
| `wsparcie_encode.png` | % retrogenów z wsparciem ENCODE per kategoria |
| `gc_polowy.png` | GC pierwsza vs druga połowa transkryptu |
| `boxplot_gc5.png` | rozkład GC pierwszej połowy między kategoriami |
| `histogram_dlugosci.png` | długość transkryptów (boxplot, skala log) |
| `gradient_gc.png` | GC wzdłuż transkryptu 5'→3' (znormalizowany) |

Wszystkie wygenerowane przez `make_charts.py`.

---

## Sekwencje i analiza GC (`strict/fasta/`)

### `intergenic_transcripts.fa` / `intronic_transcripts.fa` / `cds_transcripts.fa`
- sekwencje wybranych transkryptów (941 / 1 038 / 398).
- `>ENSG[...]|retrogen|typ|chr:start-koniec(nić)`.
- `bedtools getfasta -split -s` na `selected_encode4_*.bed`
  (eksony poskładane, nić uwzględniona).

### `gc_analysis_halves.py`
- liczy GC pierwszej i drugiej połowy każdej sekwencji
  (podział na połowy, jak Mordstein 2020) + test Wilcoxona dla par.
- czyta pliki `*_transcripts.fa`, zapisuje `gc_halves_results.csv`.

### `gc_statistics.py`
- test Shapiro-Wilk (normalność) + Mann-Whitney U (między kategoriami)
  z korektą Bonferroniego.
- czyta `gc_halves_results.csv`.

### `gc_halves_results.csv`
-  surowe dane GC.
-  kategoria, transkrypt, GC 1. połowa, GC 2. połowa, długość.
-  wygenerowany przez `gc_analysis_halves.py`.

### `selected_encode4_intergenic.bed` / `_intronic.bed` / `_cds.bed`
- wybrane transkrypty ENCODE (BED12) — kopia plików z `strict/`.

---

## Listy klasyfikacji (`strict/`)

| Plik | Co zawiera |
|------|-----------|
| `list1_intergenic_strict.bed` | 5 205 międzygenowych |
| `list2_intronic_strict.bed` | 3 191 intronowych |
| `list3_cds_one_exon_strict.bed` | 768 CDS w 1 egzonie |
| `list4_other_strict.bed` | 5 710 reszta (dotyka egzonu) |
| `list5_selfhit.bed` | 3 256 self-hit (podzbiór międzygenowych) |
| `unclassified_7.bed` | 7 niesklasyfikowanych |

Wszystkie wygenerowane przez `classify_4lists_strict.sh`.

### `single_cds.bed`
- transkrypty z dokładnie jednym egzonem kodującym (z ncbiRefSeq).
- `CDS_one_exon.py`.

### `intronic_host_transcript.bed`
- transkrypty gospodarzy dla retrogenów intronowych (pomocnicze przy klasyfikacji).
- `classify_4lists_strict.sh`.

### `selected_encode4_intergenic.bed` / `_intronic.bed` / `_cds.bed`
- wybrane transkrypty ENCODE (BED12 + nazwa retrogena w kol. 13).
- `select_encode4_transcripts.py` (+ `fix_intronic_by_ensg.py`).

### `encode4_vs_cat_list1_intergenic_strict.tsv` / `_2_intronic` / `_3_cds` / `_4_other`
- surowe wyniki przecięcia retrogenów z odczytami ENCODE.
- `validate_encode4_strict.sh` (bedtools intersect).

### `human_retrocopies_merged.bed`
-  kopia pliku wejściowego (retrokopie) w katalogu roboczym `strict/`.
-  skopiowany z głównego poziomu.

---

## Pliki NIE zawarte w ZIP (dane wejściowe, zbyt duże)

| Plik | Rola | Dlaczego nie ma |
|------|------|-----------------|
| `ncbiRefSeq.gtf` | adnotacje genów/transkryptów (hs1) | ~GB |
| `hs1.fa.gz` | sekwencja genomu hs1 | ~GB |
| `encode4_long_liftover.bed` | długie odczyty ENCODE (hg38→hs1) | ~32 MB |

Są one potrzebne do **odtworzenia** analizy, ale nie są wynikami — stąd ich brak.
