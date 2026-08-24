#!/usr/bin/env python3
"""
Analiza GC content na końcach 5' transkryptów.
Dla każdego transkryptu liczy:
  - GC w pierwszym oknie 100 nt (koniec 5')
  - GC w kolejnych oknach 100 nt (gradient 5'->3')
  - GC całego transkryptu
Porównuje 3 kategorie: intergenic, intronic, CDS.
"""

import sys
import statistics
from collections import defaultdict

WINDOW = 100  # nt
N_WINDOWS = 10  # ile okien od 5' końca

def read_fasta(path):
    seqs = {}
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(seq).upper()
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name is not None:
            seqs[name] = ''.join(seq).upper()
    return seqs

def gc_content(seq):
    if not seq:
        return 0.0
    gc = sum(1 for c in seq if c in 'GC')
    return gc / len(seq)

def window_gc(seq, win=WINDOW, n=N_WINDOWS):
    """GC content w kolejnych oknach od 5' końca."""
    result = []
    for i in range(n):
        start = i * win
        end = start + win
        if end > len(seq):
            # ostatnie okno może być krótsze
            result.append(gc_content(seq[start:]))
            break
        result.append(gc_content(seq[start:end]))
    return result

def analyze(fa_path, label):
    seqs = read_fasta(fa_path)
    data = {
        'label': label,
        'n': len(seqs),
        'gc_5prime': [],      # GC w pierwszym oknie
        'gc_overall': [],     # GC całego transkryptu
        'gc_windows': defaultdict(list),  # i-ty okienko -> lista GC
        'lengths': [],
    }
    for name, seq in seqs.items():
        w = window_gc(seq)
        data['gc_5prime'].append(w[0] if w else 0)
        data['gc_overall'].append(gc_content(seq))
        data['lengths'].append(len(seq))
        for i, g in enumerate(w):
            data['gc_windows'][i].append(g)
    return data

def mean(x):
    return sum(x) / len(x) if x else 0

def median(x):
    return statistics.median(x) if x else 0

if __name__ == '__main__':
    results = [
        analyze('intergenic_transcripts.fa', 'INTERGENIC'),
        analyze('intronic_transcripts.fa', 'INTRONIC'),
        analyze('cds_transcripts.fa', 'CDS'),
    ]

    title_5prime = "GC CONTENT NA KOŃCU 5' (okno 100 nt)"
    title_gradient = "GRADIENT GC WZDLUŻ 5'->3' (średnia w oknach 100 nt)"
    print('=' * 60)
    print('  ' + title_5prime)
    print('=' * 60)
    print()
    print(f'{"Kategoria":<15} {"n":>5} {"GC5p sr":>10} {"GC5p med":>10} {"GC calosc sr":>13}')
    print('-' * 60)
    for r in results:
        print(f'{r["label"]:<15} {r["n"]:>5} {mean(r["gc_5prime"]):>10.3f} {median(r["gc_5prime"]):>10.3f} {mean(r["gc_overall"]):>13.3f}')

    print()
    print('=' * 60)
    print('  ' + title_gradient)
    print('=' * 60)
    header = '  Okno:     ' + ''.join(f'{i*WINDOW+1:>7}' for i in range(N_WINDOWS))
    print(header)
    for r in results:
        row = f'  {r["label"]:<10}'
        for i in range(N_WINDOWS):
            g = r['gc_windows'].get(i, [])
            row += f'{mean(g):>8.3f}'
        print(row)

    # Zapisz surowe dane do CSV dla późniejszej analizy statystycznej
    import csv
    with open('gc_content_results.csv', 'w', newline='') as csvf:
        w = csv.writer(csvf)
        w.writerow(['kategoria', 'transkrypt', 'gc_5prime', 'gc_overall', 'dlugosc'])
        for r in results:
            names = list(read_fasta(f'{r["label"].lower()}_transcripts.fa').keys())
            for i, name in enumerate(names):
                w.writerow([r['label'], name, r['gc_5prime'][i], r['gc_overall'][i], r['lengths'][i]])
    print()
    print('Surowe dane zapisane do: gc_content_results.csv')
