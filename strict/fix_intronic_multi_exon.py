#!/usr/bin/env python3
"""
Dla intronowych retrogenów z wybranym transkryptem wielo-egzonowym,
spróbuj znaleźć krótki single-exon transkrypt ENCODE przez dopasowanie po nazwie genu.
Używa GTF-ów do mapowania gene_name → gene_id, potem szuka w encode4.
"""

import re
import sys
from collections import defaultdict

# 1. Zbuduj mape gene_name → gene_id z GTF-ów
print("Budowanie mapy gene_name -> gene_id...")
gene_to_ids = defaultdict(set)
for gtf_path in ['ncbiRefSeq.gtf', '../cat_liftoff/catLiftOffGenesV1.gtf']:
    try:
        with open(gtf_path) as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split('\t')
                if parts[2] != 'transcript':
                    continue
                attr = parts[8]
                gene_name = ''
                gene_id = ''
                for item in attr.split(';'):
                    item = item.strip()
                    if item.startswith('gene_name '):
                        gene_name = item.split('"')[1] if '"' in item else ''
                    elif item.startswith('gene_id '):
                        gene_id = item.split('"')[1] if '"' in item else ''
                if gene_name and gene_id:
                    gene_to_ids[gene_name].add(gene_id)
    except FileNotFoundError:
        pass
print(f"  gene_name → gene_id: {len(gene_to_ids)}")

# 2. Zbuduj mape encode4 single-exon ENSG → linia
print("Budowanie mapy encode4 single-exon...")
encode4_single = {}
with open('encode4_long_liftover.bed') as f:
    for line in f:
        if not line.strip():
            continue
        p = line.strip().split('\t')
        if int(p[9]) != 1:
            continue
        en_name = p[3]
        ensg = en_name.split('[')[0] if '[' in en_name else en_name
        if ensg not in encode4_single:
            encode4_single[ensg] = line
print(f"  single-exon ENSG: {len(encode4_single)}")

# 3. Przetworz selected_encode4_intronic.bed
print("\nPrzetwarzanie selected_encode4_intronic.bed...")
fixed = 0
kept = 0
output_lines = []

with open('selected_encode4_intronic.bed') as f:
    for line in f:
        if not line.strip():
            continue
        p = line.strip().split('\t')
        block_count = int(p[9])
        retro = p[12] if len(p) > 12 else ''

        if block_count == 1:
            output_lines.append(line)
            kept += 1
            continue

        # Multi-exon — spróbuj znaleźć single-exon przez gene_name
        parts_name = retro.split('|')
        if parts_name[0].startswith('retro_human_'):
            gn = parts_name[1] if len(parts_name) >= 2 else ''
        else:
            gn = parts_name[0]
        gn = re.sub(r'-[0-9]+$', '', gn)

        replaced = False
        if gn and gn in gene_to_ids:
            for gid in gene_to_ids[gn]:
                if gid in encode4_single:
                    new_line = encode4_single[gid].rstrip('\n')
                    new_line = new_line + '\t' + retro + '\n'
                    output_lines.append(new_line)
                    fixed += 1
                    replaced = True
                    break

        if not replaced:
            output_lines.append(line)
            kept += 1

# Zapisz
with open('selected_encode4_intronic.bed', 'w') as f:
    f.writelines(output_lines)

# Raport
multi_before = kept + fixed - sum(1 for l in output_lines if l.strip().split('\t')[9] == '1' for l in [l] )
print(f"\n=== RAPORT ===")
print(f"  Naprawione (multi-exon → single-exon przez gene_name): {fixed}")
print(f"  Pozostały bez zmian: {kept}")
print(f"  Razem: {len(output_lines)}")

# Stat
single = sum(1 for l in output_lines if int(l.strip().split('\t')[9]) == 1)
multi = len(output_lines) - single
print(f"\n  Single-exon: {single}")
print(f"  Multi-exon:  {multi}")
