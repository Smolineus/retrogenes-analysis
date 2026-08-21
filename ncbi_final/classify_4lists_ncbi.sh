#!/bin/bash
# Klasyfikacja retrogenów na 4 listy z priorytetem:
#   Lista 1: intergenowe (brak overlapu z transkryptem)
#   Lista 2: intronowe (w transkrypcie, NIE w eksonie)
#   Lista 3: CDS w 1 egzonie
#   Lista 4: eksonowe (geny wielo-egzonowe)
#
# Priorytet: CDS 1-exon > intron > ekson multi-gen > intergen
# Bez filtra self-overlapu - sprawdzamy kontekst genomowy WZGLEDEM ncbiRefSeq

set -euo pipefail

INPUT="human_retrocopies_merged.bed"
GTF="ncbiRefSeq.gtf"
MIN_OVERLAP="0.5"

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo "=== PRZYGOTOWANIE PLIKOW POMOCNICZYCH ==="

echo "[1/4] Koordynaty transkryptow (z gene_id)..."
cat "$GTF" | awk '$3=="transcript" {
    gene = $10; gsub(/"/, "", gene); gsub(/;.*/, "", gene)
    print $1"\t"($4-1)"\t"$5"\t"gene
}' > "$TMPDIR/transcripts_gene.bed"

echo "[2/4] Koordynaty transkryptow (bez gene_id)..."
cut -f1-3 "$TMPDIR/transcripts_gene.bed" > "$TMPDIR/transcripts.bed"

echo "[3/4] Koordynaty eksonow (z gene_id)..."
cat "$GTF" | awk '$3=="exon" {
    gene = $10; gsub(/"/, "", gene); gsub(/;.*/, "", gene)
    print $1"\t"($4-1)"\t"$5"\t"gene
}' > "$TMPDIR/exons_gene.bed"

echo "[4/4] Koordynaty eksonow (bez gene_id)..."
cut -f1-3 "$TMPDIR/exons_gene.bed" > "$TMPDIR/exons.bed"


echo ""
echo "=== KLASYFIKACJA ==="

extract_names() { cut -f4 "$1" | sort -u; }
sort_bed() { sort -k1,1V -k2,2n -k3,3n -k4,4 -u; }
> "$TMPDIR/assigned_names.txt"

add_gene_ids() {
    local retro_bed="$1"
    local annot_bed="$2"
    local output="$3"
    local overlap="${4:-$MIN_OVERLAP}"

    if [ ! -s "$retro_bed" ]; then
        > "$output"
        return
    fi

    bedtools intersect -a "$retro_bed" -b "$annot_bed" -wa -wb -f "$overlap" \
        | awk '
        {
            gene_id = $(NF)
            key = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6
            genekey = key "\t" gene_id
            if (!(genekey in gseen)) {
                gseen[genekey] = 1
                if (!(key in seen)) { keys[++n] = key; seen[key] = 1 }
                genes[key] = (genes[key] == "" ? gene_id : genes[key] "," gene_id)
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                print keys[i] "\t" genes[keys[i]]
            }
        }' > "$output"
}


# --- LISTA 3: CDS w 1 egzonie ---
echo "[1/4] Lista 3: CDS 1 egzon..."
bedtools intersect -a "$INPUT" -b single_cds.bed -wa -u -f "$MIN_OVERLAP" | sort_bed > "$TMPDIR/list3_classified.bed"
extract_names "$TMPDIR/list3_classified.bed" > "$TMPDIR/assigned_names.txt"

if [ -s "$TMPDIR/list3_classified.bed" ]; then
    bedtools intersect -a "$TMPDIR/list3_classified.bed" -b single_cds.bed -wa -wb -f "$MIN_OVERLAP" \
        | awk '
        {
            full_id = ""
            for (i = 10; i < NF; i++) {
                full_id = full_id (full_id ? " " : "") $i
            }
            gsub(/^gene_id: /, "", full_id)
            gsub(/;.*$/, "", full_id)
            gene = full_id

            key = $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6
            genekey = key "\t" gene
            if (!(genekey in gseen)) {
                gseen[genekey] = 1
                if (!(key in seen)) { keys[++n] = key; seen[key] = 1 }
                genes[key] = (genes[key] == "" ? gene : genes[key] "," gene)
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                print keys[i] "\t" genes[keys[i]]
            }
        }' > "$TMPDIR/list3_with_gene.bed"
else
    > "$TMPDIR/list3_with_gene.bed"
fi
echo "  -> $(cut -f4 "$TMPDIR/list3_with_gene.bed" | sort -u | wc -l) retrogenow"


# --- LISTA 2: intronowe (WYZSZY PRIORYTET niz eksonowe) ---
echo "[2/4] Lista 2: intronowe..."
# Filtr self-overlapu TYLKO dla intronowych: usun retrogeny z GTF
# zeby sprawdzic czy retrogen jest w intronie INNEGO genu
cut -f1-3 "$INPUT" > "$TMPDIR/retrocopy_mask.bed"
bedtools intersect -a "$TMPDIR/transcripts.bed" -b "$TMPDIR/retrocopy_mask.bed" -v -f 0.8 -F 0.8 > "$TMPDIR/transcripts_noself.bed"
bedtools intersect -a "$TMPDIR/exons.bed" -b "$TMPDIR/retrocopy_mask.bed" -v -f 0.8 -F 0.8 > "$TMPDIR/exons_noself.bed"
bedtools intersect -a "$TMPDIR/transcripts_gene.bed" -b "$TMPDIR/retrocopy_mask.bed" -v -f 0.8 -F 0.8 > "$TMPDIR/transcripts_gene_noself.bed"

bedtools intersect -a "$INPUT" -b "$TMPDIR/transcripts_noself.bed" -wa -u -f "$MIN_OVERLAP" | sort_bed > "$TMPDIR/tr_hits.bed"
bedtools intersect -a "$TMPDIR/tr_hits.bed" -b "$TMPDIR/exons_noself.bed" -v -f "$MIN_OVERLAP" | sort_bed > "$TMPDIR/intron_candidates.bed"
awk 'NR==FNR {assigned[$1]; next} !($4 in assigned)' "$TMPDIR/assigned_names.txt" "$TMPDIR/intron_candidates.bed" > "$TMPDIR/list2_classified.bed"
extract_names "$TMPDIR/list2_classified.bed" >> "$TMPDIR/assigned_names.txt"

add_gene_ids "$TMPDIR/list2_classified.bed" "$TMPDIR/transcripts_gene_noself.bed" "$TMPDIR/list2_with_gene.bed"
echo "  -> $(cut -f4 "$TMPDIR/list2_with_gene.bed" | sort -u | wc -l) retrogenow"


# --- LISTA 4: eksonowe ---
echo "[3/4] Lista 4: eksonowe (multi-gen)..."
bedtools intersect -a "$INPUT" -b "$TMPDIR/exons.bed" -wa -u -f "$MIN_OVERLAP" | sort_bed > "$TMPDIR/exon_hits.bed"
awk 'NR==FNR {assigned[$1]; next} !($4 in assigned)' "$TMPDIR/assigned_names.txt" "$TMPDIR/exon_hits.bed" > "$TMPDIR/list4_classified.bed"
extract_names "$TMPDIR/list4_classified.bed" >> "$TMPDIR/assigned_names.txt"

add_gene_ids "$TMPDIR/list4_classified.bed" "$TMPDIR/exons_gene.bed" "$TMPDIR/list4_with_gene.bed" "$MIN_OVERLAP"
echo "  -> $(cut -f4 "$TMPDIR/list4_with_gene.bed" | sort -u | wc -l) retrogenow"


# --- LISTA 1: intergenowe ---
echo "[4/4] Lista 1: intergenowe..."
awk 'NR==FNR {assigned[$1]; next} !($4 in assigned)' "$TMPDIR/assigned_names.txt" "$INPUT" | sort_bed > "$TMPDIR/list1_classified.bed"
awk 'BEGIN{OFS="\t"} {print $0, "intergenic"}' "$TMPDIR/list1_classified.bed" > "$TMPDIR/list1_with_gene.bed"
echo "  -> $(cut -f4 "$TMPDIR/list1_with_gene.bed" | sort -u | wc -l) retrogenow"


cp "$TMPDIR/list1_with_gene.bed" list1_intergenic.bed
cp "$TMPDIR/list2_with_gene.bed" list2_intronic.bed
cp "$TMPDIR/list3_with_gene.bed" list3_cds_one_exon.bed
cp "$TMPDIR/list4_with_gene.bed" list4_exonic_multigene.bed


echo ""
echo "========================================"
echo "            PODSUMOWANIE"
echo "========================================"
printf "%-35s %8s\n" "Plik" "Unikalne"

declare -A totals
for f in list1_intergenic.bed list2_intronic.bed list3_cds_one_exon.bed list4_exonic_multigene.bed; do
    count=$(cut -f4 "$f" 2>/dev/null | sort -u | wc -l)
    printf "%-35s %8d\n" "$f" "$count"
    totals["$f"]=$count
done

sum=$(( totals[list1_intergenic.bed] + totals[list2_intronic.bed] + totals[list3_cds_one_exon.bed] + totals[list4_exonic_multigene.bed] ))
echo "----------------------------------------"
printf "%-35s %8d\n" "SUMA (4 listy)" "$sum"
printf "%-35s %8d\n" "Oryginalny plik" "$(wc -l < "$INPUT")"

echo ""
echo "Pliki wynikowe:"
echo "  list1_intergenic.bed        - międzygenowe"
echo "  list2_intronic.bed          - intronowe"
echo "  list3_cds_one_exon.bed      - CDS w 1 egzonie"
echo "  list4_exonic_multigene.bed  - eksonowe (geny wielo-egzonowe)"
echo ""
