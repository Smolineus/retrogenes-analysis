#!/usr/bin/env python3
"""
Selekcja transkryptów ENCODE dla retrogenów — jeden transkrypt na retrogen.
Priorytet: 1) single-exon  2) długość najbliższa retrogenowi  3) lider 5'.
Dla listy intronowej: dodatkowy filtr długości transkryptu ≤ 2 × długość retrogenu.
"""

import sys
from collections import defaultdict

def parse_tsv(tsv_path, retrogene_list_name, max_length_ratio=None):
    retro_groups = defaultdict(list)
    total_matches = 0
    retrogenes_with_matches = set()

    with open(tsv_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 20:
                continue

            retro_start  = int(parts[13])
            retro_end    = int(parts[14])
            retro_name   = parts[15]

            retrogene_len = retro_end - retro_start

            en_strand = parts[5]
            en_start  = int(parts[1])
            en_end    = int(parts[2])

            # Transcript length = sum of blockSizes
            block_sizes_str = parts[10].rstrip(',').split(',')
            transcript_len = sum(int(x) for x in block_sizes_str if x)

            # Length filter (for intronic)
            if max_length_ratio is not None:
                if transcript_len > max_length_ratio * retrogene_len:
                    continue

            # Single-exon check
            block_count = int(parts[9])
            is_single_exon = (block_count == 1)

            # 5' leader
            if en_strand == '+':
                leader = retro_start - en_start
            else:
                leader = en_end - retro_end
            if leader < 0:
                leader = 0

            # Distance to retrogene length (smaller = better)
            distance = abs(transcript_len - retrogene_len)

            # Store: (parts, is_single_exon, distance, leader, transcript_len)
            retro_groups[retro_name].append((parts, is_single_exon, distance, leader, transcript_len))
            retrogenes_with_matches.add(retro_name)
            total_matches += 1

    # Select best transcript per retrogene
    selected = []
    for name, candidates in retro_groups.items():
        if max_length_ratio is not None:
            # Intronic: single-exon ma pierwszeństwo (odróżnia retrogen od gospodarza)
            # sort: single_exon DESC, distance ASC, leader DESC
            candidates.sort(key=lambda x: (x[1], -x[2], x[3]), reverse=True)
        else:
            # Intergenic/CDS: brak preferencji single-exon
            # sort: distance ASC, leader DESC
            candidates.sort(key=lambda x: (-x[2], x[3]), reverse=True)
        best = candidates[0]
        encode4_parts = best[0][:12]
        selected.append((encode4_parts, name, best[3], best[4]))

    # Report
    print(f"=== {retrogene_list_name} ===")
    print(f"  Total TSV lines:          {total_matches}")
    print(f"  Unique retrogenes:        {len(retrogenes_with_matches)}")
    print(f"  Selected (output):        {len(selected)}")
    if max_length_ratio:
        filtered_out = sum(1 for name in retrogenes_with_matches if name not in [s[1] for s in selected])
        print(f"  Filtered out (length):    {filtered_out}")
    print()

    return selected


def write_bed12(selected, output_path):
    with open(output_path, 'w') as out:
        for encode4_parts, retro_name, leader, trans_len in selected:
            bed12 = encode4_parts[:12]
            line = '\t'.join(bed12) + '\t' + retro_name
            out.write(line + '\n')


def process_list(tsv_path, output_path, label, max_length_ratio=None):
    selected = parse_tsv(tsv_path, label, max_length_ratio)
    write_bed12(selected, output_path)
    return len(selected)


if __name__ == '__main__':
    n1 = process_list(
        'encode4_vs_cat_list1_intergenic_strict.tsv',
        'selected_encode4_intergenic.bed',
        'INTERGENIC'
    )

    n2 = process_list(
        'encode4_vs_cat_list2_intronic_strict.tsv',
        'selected_encode4_intronic.bed',
        'INTRONIC',
        max_length_ratio=2.0
    )

    n3 = process_list(
        'encode4_vs_cat_list3_cds_one_exon_strict.tsv',
        'selected_encode4_cds.bed',
        'CDS 1 EXON'
    )

    print("========================================")
    print("            PODSUMOWANIE")
    print("========================================")
    print(f"  Intergenic: {n1} transkryptów")
    print(f"  Intronic:   {n2} transkryptów")
    print(f"  CDS 1 egzon:{n3} transkryptów")
    print(f"  RAZEM:      {n1 + n2 + n3}")
    print()
    print("Pliki wyjściowe (BED12 + nazwa retrogenu w kol. 13):")
    print("  selected_encode4_intergenic.bed")
    print("  selected_encode4_intronic.bed")
    print("  selected_encode4_cds.bed")
