#!/bin/bash
# Klasyfikacja retrogenów na 4 listy z priorytetem:
#   Lista 1: intergenowe
#   Lista 2: intronowe
#   Lista 3: CDS w 1 egzonie
#   Lista 4: eksonowe (geny wielo-egzonowe)
#
# Priorytet: CDS 1-exon > ekson multi-gen > intron > intergen
# Prog overlapu: 50% retrogenu (-f 0.5) dla checkow pozytywnych
#
# Wyjscie: BED6 + gene_id (kolumna 7), wiele genow laczone przecinkiem

set -euo pipefail

INPUT="human_retrocopies_merged.bed"
GTF="hs1.ncbiRefSeq.bigZip.gtf.gz"
MIN_OVERLAP="0.5"

TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo "=== PRZYGOTOWANIE PLIKOW POMOCNICZYCH ==="

echo "[1/4] Koordynaty transkryptow (z gene_id)..."
zcat "$GTF" | awk '$3=="transcript" {
    gene = $10; gsub(/"/, "", gene); gsub(/;.*/, "", gene)
    print $1"\t"($4-1)"\t"$5"\t"gene
}' > "$TMPDIR/transcripts_gene.bed"

echo "[2/4] Koordynaty transkryptow (bez gene_id, do klasyfikacji)..."
cut -f1-3 "$TMPDIR/transcripts_gene.bed" > "$TMPDIR/transcripts.bed"

echo "[3/4] Koordynaty eksonow (z gene_id)..."
zcat "$GTF" | awk '$3=="exon" {
    gene = $10; gsub(/"/, "", gene); gsub(/;.*/, "", gene)
    print $1"\t"($4-1)"\t"$5"\t"gene
}' > "$TMPDIR/exons_gene.bed"

echo "[4/4] Koordynaty eksonow (bez gene_id, do klasyfikacji)..."
cut -f1-3 "$TMPDIR/exons_gene.bed" > "$TMPDIR/exons.bed"

if [ ! -f single_cds.bed ]; then
    echo "Uruchamiam CDS_one_exon.py..."
    python3 CDS_one_exon.py
fi


echo ""
echo "=== KLASYFIKACJA ==="

extract_names() { cut -f4 "$1" | sort -u; }
# sort_bed: sortuje BED wg chr (wersyjnie), start, end, name
sort_bed() { sort -k1,1V -k2,2n -k3,3n -k4,4 -u; }
> "$TMPDIR/assigned_names.txt"

# Funkcja: dodaje gene_id z pliku adnotacji do retrogenow z listy
# $1 = plik z retrogenami (BED6)
# $2 = plik z adnotacjami (BED4: chr,start,end,gene_id)
# $3 = plik wyjsciowy (BED7: chr,start,end,name,score,strand,gene_ids)
add_gene_ids() {
    local retro_bed="$1"
    local annot_bed="$2"
    local output="$3"

    if [ ! -s "$retro_bed" ]; then
        > "$output"
        return
    fi

    bedtools intersect -a "$retro_bed" -b "$annot_bed" -wa -wb -f "$MIN_OVERLAP" \
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

# --- LISTA 3: CDS w 1 egzonie (najwyzszy priorytet) ---
echo "[1/4] Lista 3: CDS 1 egzon..."
bedtools intersect -a "$INPUT" -b single_cds.bed -wa -u -f "$MIN_OVERLAP" | sort_bed > "$TMPDIR/list3_classified.bed"
extract_names "$TMPDIR/list3_classified.bed" > "$TMPDIR/assigned_names.txt"

# Dodajemy gene_id z single_cds (kolumna 4 ma format: gene_id: X; transcript_id: Y)
if [ -s "$TMPDIR/list3_classified.bed" ]; then
    bedtools intersect -a "$TMPDIR/list3_classified.bed" -b single_cds.bed -wa -wb -f "$MIN_OVERLAP" \
        | awk '
        {
            # Rekonstruujemy full_id z pol $10..$(NF-1) (rozbite spacjami)
            full_id = ""
            for (i = 10; i < NF; i++) {
                full_id = full_id (full_id ? " " : "") $i
            }
            # Wyciagamy gene_id: "gene_id: SPRY3; transcript_id: NM_001304990.2" -> "SPRY3"
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


# --- LISTA 4: eksonowe genow wielo-egzonowych ---
echo "[2/4] Lista 4: eksonowe (multi-gen)..."
bedtools intersect -a "$INPUT" -b "$TMPDIR/exons.bed" -wa -u -f "$MIN_OVERLAP" | sort_bed > "$TMPDIR/exon_hits.bed"
awk 'NR==FNR {assigned[$1]; next} !($4 in assigned)' "$TMPDIR/assigned_names.txt" "$TMPDIR/exon_hits.bed" > "$TMPDIR/list4_classified.bed"
extract_names "$TMPDIR/list4_classified.bed" >> "$TMPDIR/assigned_names.txt"

add_gene_ids "$TMPDIR/list4_classified.bed" "$TMPDIR/exons_gene.bed" "$TMPDIR/list4_with_gene.bed"
echo "  -> $(cut -f4 "$TMPDIR/list4_with_gene.bed" | sort -u | wc -l) retrogenow"


# --- LISTA 2: intronowe (w transkrypcie, NIE w eksonie) ---
echo "[3/4] Lista 2: intronowe..."
bedtools intersect -a "$INPUT" -b "$TMPDIR/transcripts.bed" -wa -u -f "$MIN_OVERLAP" | sort_bed > "$TMPDIR/tr_hits.bed"
bedtools intersect -a "$TMPDIR/tr_hits.bed" -b "$TMPDIR/exons.bed" -v | sort_bed > "$TMPDIR/intron_candidates.bed"
awk 'NR==FNR {assigned[$1]; next} !($4 in assigned)' "$TMPDIR/assigned_names.txt" "$TMPDIR/intron_candidates.bed" > "$TMPDIR/list2_classified.bed"
extract_names "$TMPDIR/list2_classified.bed" >> "$TMPDIR/assigned_names.txt"

add_gene_ids "$TMPDIR/list2_classified.bed" "$TMPDIR/transcripts_gene.bed" "$TMPDIR/list2_with_gene.bed"
echo "  -> $(cut -f4 "$TMPDIR/list2_with_gene.bed" | sort -u | wc -l) retrogenow"


# --- LISTA 1: intergenowe (reszta) ---
echo "[4/4] Lista 1: intergenowe..."
awk 'NR==FNR {assigned[$1]; next} !($4 in assigned)' "$TMPDIR/assigned_names.txt" "$INPUT" | sort_bed > "$TMPDIR/list1_classified.bed"
awk 'BEGIN{OFS="\t"} {print $0, "intergenic"}' "$TMPDIR/list1_classified.bed" > "$TMPDIR/list1_with_gene.bed"
echo "  -> $(cut -f4 "$TMPDIR/list1_with_gene.bed" | sort -u | wc -l) retrogenow"


# --- KOPIUJ WYNIKI ---
cp "$TMPDIR/list1_with_gene.bed" list1_intergenic.bed
cp "$TMPDIR/list2_with_gene.bed" list2_intronic.bed
cp "$TMPDIR/list3_with_gene.bed" list3_cds_one_exon.bed
cp "$TMPDIR/list4_with_gene.bed" list4_exonic_multigene.bed


# --- RAPORT ---
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
echo "Format plikow wynikowych (BED7):"
echo "  chr start end retro_name score strand gene_id(gene_id,...)"
echo ""
echo "Pliki wynikowe:"
echo "  list1_intergenic.bed        - międzygenowe"
echo "  list2_intronic.bed          - intronowe"
echo "  list3_cds_one_exon.bed      - CDS w 1 egzonie"
echo "  list4_exonic_multigene.bed  - eksonowe (geny wielo-egzonowe)"
echo ""
