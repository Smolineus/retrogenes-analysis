#!/usr/bin/env python3
"""
Dla intronowych retrogenow z wybranym transkryptem multi-exon,
znajdz single-exon transkrypt z tego samego ENSG nachodzacy w JAKIKOLWIEK stopniu.
"""
import sys
from collections import defaultdict

# 1. Zbuduj mape ENSG → lista single-exon entries z encode4
print("Budowanie mapy encode4 single-exon...")
encode4_single = defaultdict(list)
with open('encode4_long_liftover.bed') as f:
    for line in f:
        if not line.strip(): continue
        p = line.strip().split('\t')
        if int(p[9]) != 1: continue
        en_name = p[3]
        ensg = en_name.split('[')[0] if '[' in en_name else en_name
        encode4_single[ensg].append(line)

print(f"  Unikalnych ENSG z single-exon: {len(encode4_single)}")

# 2. Wczytaj selected intronic, znajdz multi-exon
print("\nPrzetwarzanie...")
fixed = 0
not_fixed = []
output_lines = []

with open('selected_encode4_intronic.bed') as f:
    for line in f:
        if not line.strip(): continue
        p = line.strip().split('\t')
        block_count = int(p[9])
        retro = p[12] if len(p) > 12 else ''

        if block_count == 1:
            output_lines.append(line)
            continue

        # Multi-exon: szukaj single-exon z tego samego ENSG nachodzacego na retrogen
        en_name = p[3]
        ensg = en_name.split('[')[0] if '[' in en_name else en_name
        
        # Znajdz koordynaty retrogenu
        retro_coords = None
        with open('human_retrocopies_merged.bed') as rfile:
            for rl in rfile:
                rp = rl.strip().split('\t')
                if rp[3] == retro:
                    retro_coords = (rp[0], int(rp[1]), int(rp[2]))
                    break

        replaced = False
        if retro_coords and ensg in encode4_single:
            for se_line in encode4_single[ensg]:
                se_p = se_line.strip().split('\t')
                se_chr = se_p[0]
                se_start = int(se_p[1])
                se_end = int(se_p[2])
                r_chr, r_start, r_end = retro_coords
                
                # Sprawdz overlap (jakikolwiek)
                if se_chr == r_chr and se_start < r_end and se_end > r_start:
                    new_line = se_line.rstrip('\n') + '\t' + retro + '\n'
                    output_lines.append(new_line)
                    fixed += 1
                    replaced = True
                    break

        if not replaced:
            output_lines.append(line)
            not_fixed.append((en_name, retro, block_count))

with open('selected_encode4_intronic.bed', 'w') as f:
    f.writelines(output_lines)

print(f"\nNaprawione: {fixed}")
print(f"Pozostałe multi-exon: {len(not_fixed)}")
if not_fixed:
    print("\nPozostałe:")
    for ensg, retro, bc in not_fixed:
        print(f"  {ensg:45s} → {retro}")

# Stat
single = sum(1 for l in output_lines if int(l.strip().split('\t')[9]) == 1)
multi = len(output_lines) - single
print(f"\nSingle-exon: {single}")
print(f"Multi-exon:  {multi}")
