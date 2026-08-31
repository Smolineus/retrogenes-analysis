#!/usr/bin/env python3
"""
Analiza GC content — podział transkryptu na połowy (jak Mordstein 2020).
Dla każdego transkryptu liczy GC pierwszej połowy (5') i drugiej połowy (3'),
a następnie testem Wilcoxona dla par (signed-rank) sprawdza, czy pierwsza
połowa ma istotnie wyższe GC.
"""

import statistics
from collections import defaultdict
from scipy import stats


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


def halves_gc(seq):
    """Dzieli sekwencję na dwie równe połowy i liczy GC każdej."""
    mid = len(seq) // 2
    first_half = seq[:mid]
    second_half = seq[mid:]
    return gc_content(first_half), gc_content(second_half)


def analyze(fa_path, label):
    seqs = read_fasta(fa_path)
    data = {
        'label': label,
        'n': len(seqs),
        'gc_first': [],   # GC pierwszej połowy (5')
        'gc_second': [],  # GC drugiej połowy (3')
        'lengths': [],
    }
    for name, seq in seqs.items():
        gc1, gc2 = halves_gc(seq)
        data['gc_first'].append(gc1)
        data['gc_second'].append(gc2)
        data['lengths'].append(len(seq))
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

    print('=' * 70)
    print("  GC CONTENT — PODZIAŁ NA POŁOWY (pierwsza vs druga połowa)")
    print('=' * 70)
    print()
    print(f'{"Kategoria":<15} {"n":>5} {"GC 1. połowa":>13} {"GC 2. połowa":>13} {"różnica":>9}')
    print('-' * 70)
    for r in results:
        diff = mean(r['gc_first']) - mean(r['gc_second'])
        print(f'{r["label"]:<15} {r["n"]:>5} '
              f'{mean(r["gc_first"]):>13.3f} {mean(r["gc_second"]):>13.3f} {diff:>+9.3f}')

    print()
    print('=' * 70)
    print('  TEST WILCOXONA DLA PAR (pierwsza vs druga połowa)')
    print('=' * 70)
    print()
    for r in results:
        w_stat, p = stats.wilcoxon(r['gc_first'], r['gc_second'], alternative='two-sided')
        # kierunek: czy pierwsza połowa jest WYŻSZA
        w_one, p_one = stats.wilcoxon(r['gc_first'], r['gc_second'], alternative='greater')
        sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'n.s.'))
        direction = '1. połowa WYŻSZA' if p_one < 0.05 else ('2. połowa WYŻSZA' if mean(r['gc_first']) < mean(r['gc_second']) else 'brak różnicy')
        print(f'  {r["label"]:<15} W={w_stat:.1f} p={p:.3e} {sig}  ({direction})')

    print()
    # Zapisz surowe dane
    import csv
    with open('gc_halves_results.csv', 'w', newline='') as csvf:
        w = csv.writer(csvf)
        w.writerow(['kategoria', 'transkrypt', 'gc_pierwsza_polowa', 'gc_druga_polowa', 'dlugosc'])
        for r in results:
            names = list(read_fasta(f'{r["label"].lower()}_transcripts.fa').keys())
            for i, name in enumerate(names):
                w.writerow([r['label'], name, r['gc_first'][i], r['gc_second'][i], r['lengths'][i]])
    print()
    print('Surowe dane zapisane do: gc_halves_results.csv')
