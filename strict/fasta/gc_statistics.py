#!/usr/bin/env python3
"""
Analiza statystyczna GC content (pierwsza połowa transkryptu, 5') między kategoriami.

Kroki:
1. Test normalności Shapiro-Wilka — uzasadnia wybór testu nieparametrycznego
   (jeśli dane NIE są normalne -> Mann-Whitney U)
2. Test Manna-Whitneya U — porównanie GC 5' między kategoriami
   (z korektą Bonferroniego dla wielokrotnych porównań)
"""

import csv
from collections import defaultdict
from scipy import stats


def load_csv(path):
    data = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data[row['kategoria']].append(float(row['gc_pierwsza_polowa']))
    return data


if __name__ == '__main__':
    data = load_csv('gc_halves_results.csv')
    cats = ['INTERGENIC', 'INTRONIC', 'CDS']

    # Liczba porównań parami (do korekty Bonferroniego)
    n_comparisons = len(cats) * (len(cats) - 1) // 2
    bonferroni_alpha = 0.05 / n_comparisons

    print('=' * 64)
    print('  1. TEST NORMALNOŚCI (Shapiro-Wilk) — GC 1. połowy (5\')')
    print('=' * 64)
    print()
    print(f'{"Kategoria":<15} {"n":>6} {"W":>10} {"p-value":>12} {"Wniosek":<20}')
    print('-' * 64)
    for c in cats:
        w_stat, p = stats.shapiro(data[c])
        concl = 'NIE normalny' if p < 0.05 else 'normalny (n.s.)'
        print(f'{c:<15} {len(data[c]):>6} {w_stat:>10.4f} {p:>12.3e} {concl:<20}')
    print()
    print('Wniosek: jeśli p < 0.05, rozkład nie jest normalny -> test nieparametryczny.')
    print()

    print('=' * 64)
    print('  2. TEST MANN-WHITNEY U — GC 1. połowy (5\')')
    print('=' * 64)
    print()
    print(f'Korekta Bonferroniego: α = 0.05/{n_comparisons} = {bonferroni_alpha:.4f}')
    print('(wynik istotny, gdy p < skorygowane α)')
    print()
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            a, b = cats[i], cats[j]
            u, p = stats.mannwhitneyu(data[a], data[b], alternative='two-sided')
            sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'n.s.'))
            # Czy istotne po korekcie Bonferroniego?
            bonf = 'TAK' if p < bonferroni_alpha else 'nie'
            print(f'  {a:<12} vs {b:<12} n1={len(data[a])} n2={len(data[b])} '
                  f'U={u:.1f} p={p:.3e}  {sig}  (po Bonf.: {bonf})')
    print()
