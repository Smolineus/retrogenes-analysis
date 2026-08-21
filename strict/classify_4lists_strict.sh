#!/bin/bash
# ŚCISLA klasyfikacja retrogenów:
#   Lista 1: 100% międzygenowe (ZERO overlapu z transkryptem)
#   Lista 2: 100% intronowe (w transkrypcie, ZERO overlapu z eksonem)
#   Lista 3: CDS w 1 egzonie (≥50% overlap z single-CDS)
#   Lista 4: reszta (nie spełniają kryteriów 1-3)
#
# Priorytet: CDS 1-exon > intron (0% ekson) > ekson multi-gen > intergen (0% transkrypt)

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
# Filtr self-overlapu PO NAZWIE GENU: usuń z GTF transkrypty/eksony
# których gene_name pasuje do nazwy retrogenu
awk -F'\t' '{print $4}' "$INPUT" | awk -F'|' '{
    if ($1 ~ /^retro_human_/) {
        if (NF >= 2) name=$2; else name=""
    } else name=$1
    gsub(/-[0-9]+$/, "", name)
    if (name != "") print name
}' | sort -u > "$TMPDIR/retrocopy_gene_names.txt"

awk -F'\t' 'NR==FNR {
    names[$1]=1; next
}
{
    gene = $4
    if (!(gene in names)) print
}' "$TMPDIR/retrocopy_gene_names.txt" "$TMPDIR/exons_gene.bed" | cut -f1-3 > "$TMPDIR/exons_noself.bed"

# Transkryptow NIE filtrujemy dla overlap-check - retrogen musi znalezc gospodarza
# nawet jesli inny retrogen dzieli z nim nazwe genu
cp "$TMPDIR/transcripts.bed" "$TMPDIR/transcripts_noself.bed"

# Ale dla gene_id lookup filtrujemy - zeby dostac nazwe gospodarza, nie siebie
awk -F'\t' 'NR==FNR {
    names[$1]=1; next
}
{
    gene = $4
    if (!(gene in names)) print
}' "$TMPDIR/retrocopy_gene_names.txt" "$TMPDIR/transcripts_gene.bed" > "$TMPDIR/transcripts_gene_noself.bed"

bedtools intersect -a "$INPUT" -b "$TMPDIR/transcripts_noself.bed" -wa -u -f "$MIN_OVERLAP" | sort_bed > "$TMPDIR/tr_hits.bed"
bedtools intersect -a "$TMPDIR/tr_hits.bed" -b "$TMPDIR/exons_noself.bed" -v | sort_bed > "$TMPDIR/intron_candidates.bed"
awk 'NR==FNR {assigned[$1]; next} !($4 in assigned)' "$TMPDIR/assigned_names.txt" "$TMPDIR/intron_candidates.bed" > "$TMPDIR/list2_classified.bed"
extract_names "$TMPDIR/list2_classified.bed" >> "$TMPDIR/assigned_names.txt"

add_gene_ids "$TMPDIR/list2_classified.bed" "$TMPDIR/transcripts_gene.bed" "$TMPDIR/list2_with_gene.bed"
echo "  -> $(cut -f4 "$TMPDIR/list2_with_gene.bed" | sort -u | wc -l) retrogenow (przed korekta)"

# Korekta: retrogeny z listy 2 ktore maja TYLKO siebie jako gene_id → lista 1
awk -F'\t' '{
    split($4, parts, "|")
    if (parts[1] ~ /^retro_human_/) n=(NF>=2 ? parts[2] : "")
    else n=parts[1]
    gsub(/-[0-9]+$/, "", n)
    split($7, genes, ",")
    has_other=0
    for(g in genes) {
        gsub(/^[ \t]+|[ \t]+$/,"",genes[g])
        if(genes[g]!=n) has_other=1
    }
    if(!has_other && n!="") print
}' "$TMPDIR/list2_with_gene.bed" > "$TMPDIR/list2_to_intergenic.bed"

awk -F'\t' '{
    split($4, parts, "|")
    if (parts[1] ~ /^retro_human_/) n=(NF>=2 ? parts[2] : "")
    else n=parts[1]
    gsub(/-[0-9]+$/, "", n)
    split($7, genes, ",")
    has_other=0
    for(g in genes) {
        gsub(/^[ \t]+|[ \t]+$/,"",genes[g])
        if(genes[g]!=n) has_other=1
    }
    if(has_other || n=="") print
}' "$TMPDIR/list2_with_gene.bed" > "$TMPDIR/list2_clean.bed"

mv "$TMPDIR/list2_clean.bed" "$TMPDIR/list2_with_gene.bed"

# Druga korekta: koordynatowa — retrogeny >=90% na ekson,
# wszystkie gene_id to nazwy retrogenow → self-hit → lista 1
bedtools intersect -a "$TMPDIR/list2_with_gene.bed" -b "$TMPDIR/exons.bed" -wa -wb -f 0.9 -F 0.9 \
    > "$TMPDIR/list2_exon_overlap.tsv"

python3 -c "
import re

# Wczytaj nazwy retrogenow
rnames=set()
with open('$TMPDIR/retrocopy_gene_names.txt') as f:
    for l in f:
        n=l.strip()
        if n: rnames.add(n)

# Sprawdz ktore retrogeny z overlapem maja same retrogeny w gene_id
to_selfhit=set()
with open('$TMPDIR/list2_exon_overlap.tsv') as f:
    for l in f:
        p=l.strip().split('\t')
        name=p[3]
        gids=set(g.strip() for g in p[6].split(',') if g.strip())
        if gids and all(g in rnames for g in gids):
            to_selfhit.add(name)

# Podziel list2
with open('$TMPDIR/list2_with_gene.bed') as fin, \
     open('$TMPDIR/list2_clean2.bed','w') as fclean, \
     open('$TMPDIR/list2_coord_selfhit.bed','w') as fself:
    for l in fin:
        name=l.strip().split('\t')[3]
        if name in to_selfhit: fself.write(l)
        else: fclean.write(l)

print(f'{len(to_selfhit)} koordynatowych self-hit')
" 2>/dev/null

mv "$TMPDIR/list2_clean2.bed" "$TMPDIR/list2_with_gene.bed" 2>/dev/null
cat "$TMPDIR/list2_coord_selfhit.bed" >> "$TMPDIR/list2_to_intergenic.bed" 2>/dev/null

cp "$TMPDIR/list2_to_intergenic.bed" list5_selfhit.bed
echo "  -> po korekcie: $(cut -f4 "$TMPDIR/list2_with_gene.bed" 2>/dev/null | sort -u | wc -l) intronowe, $(wc -l < "$TMPDIR/list2_to_intergenic.bed" 2>/dev/null) self-hit"


# --- LISTA 4: reszta (NIE-spelniajace 1-3) ---
echo "[3/4] Lista 4: reszta (eksonowe / graniczne)..."
# Bez progu -f: kazdy overlap z eksonem (nawet 1bp) → lista 4
# zeby nic nie wpadlo w luke miedzy -v listy 2 a -f 0.1
bedtools intersect -a "$INPUT" -b "$TMPDIR/exons.bed" -wa -u | sort_bed > "$TMPDIR/exon_hits.bed"
awk 'NR==FNR {assigned[$1]; next} !($4 in assigned)' "$TMPDIR/assigned_names.txt" "$TMPDIR/exon_hits.bed" > "$TMPDIR/list4_classified.bed"
extract_names "$TMPDIR/list4_classified.bed" >> "$TMPDIR/assigned_names.txt"

bedtools intersect -a "$TMPDIR/list4_classified.bed" -b "$TMPDIR/exons_gene.bed" -wa -wb | awk '{
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
}' > "$TMPDIR/list4_with_gene.bed"
echo "  -> $(cut -f4 "$TMPDIR/list4_with_gene.bed" | sort -u | wc -l) retrogenow"


# --- LISTA 1: intergenowe ---
echo "[4/4] Lista 1: intergenowe..."
awk 'NR==FNR {assigned[$1]; next} !($4 in assigned)' "$TMPDIR/assigned_names.txt" "$INPUT" | sort_bed > "$TMPDIR/list1_classified.bed"
awk 'BEGIN{OFS="\t"} {print $0, "intergenic"}' "$TMPDIR/list1_classified.bed" > "$TMPDIR/list1_with_gene.bed"
# Dolacz retrogeny przeniesione z listy 2 (self-only)
cat "$TMPDIR/list2_to_intergenic.bed" >> "$TMPDIR/list1_with_gene.bed"
echo "  -> $(cut -f4 "$TMPDIR/list1_with_gene.bed" | sort -u | wc -l) retrogenow"


cp "$TMPDIR/list1_with_gene.bed" list1_intergenic_strict.bed
cp "$TMPDIR/list2_with_gene.bed" list2_intronic_strict.bed
cp "$TMPDIR/list3_with_gene.bed" list3_cds_one_exon_strict.bed
cp "$TMPDIR/list4_with_gene.bed" list4_other_strict.bed


echo ""
echo "========================================"
echo "            PODSUMOWANIE"
echo "========================================"
printf "%-35s %8s\n" "Plik" "Unikalne"

declare -A totals
for f in list1_intergenic_strict.bed list2_intronic_strict.bed list3_cds_one_exon_strict.bed list4_other_strict.bed; do
    count=$(cut -f4 "$f" 2>/dev/null | sort -u | wc -l)
    printf "%-35s %8d\n" "$f" "$count"
    totals["$f"]=$count
done

sum=$(( totals[list1_intergenic_strict.bed] + totals[list2_intronic_strict.bed] + totals[list3_cds_one_exon_strict.bed] + totals[list4_other_strict.bed] ))
echo "----------------------------------------"
printf "%-35s %8d\n" "SUMA (4 listy)" "$sum"
printf "%-35s %8d\n" "Oryginalny plik" "$(wc -l < "$INPUT")"

echo ""
echo "Pliki wynikowe:"
echo "  list1_intergenic_strict.bed        - międzygenowe"
echo "  list2_intronic_strict.bed          - intronowe"
echo "  list3_cds_one_exon_strict.bed      - CDS w 1 egzonie"
echo "  list4_other_strict.bed  - eksonowe (geny wielo-egzonowe)"
echo ""
