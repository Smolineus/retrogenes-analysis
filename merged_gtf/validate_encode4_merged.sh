#!/bin/bash
# Walidacja 4 list retrogenow z dlugimi odczytami ENCODE (po liftoverze na hs1)
# Sprawdza, ktore retrogeny pokrywaja sie z encode4_long_liftover.bed

set -euo pipefail

ENCODE="encode4_long_liftover.bed"
MIN_F_A="0.1"
MIN_F_B="0.1"

echo "========================================"
echo "  WALIDACJA ENCODE4 vs 4 LISTY"
echo "========================================"
echo ""

declare -A LISTS
LISTS=(
    ["list1_intergenic_merged.bed"]="Intergenowe"
    ["list2_intronic_merged.bed"]="Intronowe"
    ["list3_cds_one_exon_merged.bed"]="CDS 1 egzon"
    ["list4_exonic_multigene_merged.bed"]="Eksonowe multi-gen"
)

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

printf "%-30s %10s %10s %10s\n" "Kategoria" "W liscie" "W ENCODE" "% wsparcia"
echo "--------------------------------------------------------------------"

for bed in "${!LISTS[@]}"; do
    label="${LISTS[$bed]}"
    tsv_out="encode4_vs_cat_${bed%.bed}.tsv"

    if [ ! -f "$bed" ]; then
        echo "BRAK PLIKU $bed"
        continue
    fi

    # Liczba unikalnych retrogenow w liscie
    total=$(cut -f4 "$bed" | sort -u | wc -l)

    # Intersect encode4 z lista
    bedtools intersect -a "$ENCODE" -b "$bed" -wo -f "$MIN_F_A" -F "$MIN_F_B" \
        > "$TMPDIR/intersect.tsv" 2>/dev/null

    # Wyciagamy unikalne nazwy retrogenow, ktore maja overlap z encode4
    # encode4_long_liftover.bed: kolumna 4 to nazwa transkryptu ENCODE
    # nasze listy: w kolumnie 4 (w pliku B) bedzie nazwa retrogenu
    # bedtools -wo dodaje kolumny B na koncu -> nazwa retrogenu jest w $(NF-...)
    # Prosciej: uzyc bedtools z odwrotna kolejnocia A/B zeby dostac nazwy retrogenow
    bedtools intersect -b "$ENCODE" -a "$bed" -wa -u -f "$MIN_F_A" -F "$MIN_F_B" \
        > "$TMPDIR/validated.bed" 2>/dev/null

    validated=$(cut -f4 "$TMPDIR/validated.bed" | sort -u | wc -l)

    # Procent wsparcia
    if [ "$total" -gt 0 ]; then
        pct=$(awk "BEGIN {printf \"%.1f\", ($validated/$total)*100}")
    else
        pct="0.0"
    fi

    printf "%-30s %10d %10d %10s%%\n" "$label" "$total" "$validated" "$pct"

    # Zachowaj plik tsv z pelnym intersectem (A=encode4, B=lista, -wo)
    bedtools intersect -a "$ENCODE" -b "$bed" -wo -f "$MIN_F_A" -F "$MIN_F_B" \
        > "$tsv_out" 2>/dev/null
    echo "  -> $tsv_out"
done

echo ""
echo "Pliki wyjsciowe (*.tsv): encode4 (A) vs lista (B), format: kolumny A + kolumny B + overlap_bp"
echo "Do wyciagniecia nazw retrogenow z tsva: cut -f\$(head -1 encode4_long_liftover.bed | awk '{print NF+4}') plik.tsv | sort -u"
echo ""
