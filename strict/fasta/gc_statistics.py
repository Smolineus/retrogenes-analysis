#!/usr/bin/env python3
"""
Analiza statystyczna GC content (koniec 5') między kategoriami.

Kroki:
1. Test normalności Andersona-Darlinga — uzasadnia wybór testu nieparametrycznego
   (jeśli dane NIE są normalne -> Mann-Whitney U)
2. Test Manna-Whitneya U — porównanie GC 5' między kategoriami
"""

import csv
import math
import statistics
from collections import defaultdict, Counter


def load_csv(path):
    data = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data[row['kategoria']].append(float(row['gc_5prime']))
    return data


def anderson_darling(x):
    """Test normalności Andersona-Darlinga.
    H0: dane pochodzą z rozkładu normalnego.
    Zwraca (A^2, p-value). p < 0.05 -> odrzucamy normalność."""
    n = len(x)
    if n < 8:
        return None, None
    mean = statistics.mean(x)
    sd = statistics.stdev(x)
    if sd == 0:
        return 0.0, 1.0

    x_sorted = sorted(x)
    z = [(xi - mean) / sd for xi in x_sorted]

    def norm_cdf(t):
        return 0.5 * (1 + math.erf(t / math.sqrt(2)))

    A2 = -n
    for i in range(1, n + 1):
        zi = z[i - 1]
        zn1i = z[n - i]
        A2 -= (2 * i - 1) * (math.log(norm_cdf(zi) + 1e-300) + math.log(1 - norm_cdf(zn1i) + 1e-300)) / n

    # Korekta dla skończonego n
    A2_star = A2 * (1 + 0.75 / n + 2.25 / (n * n))

    # Przybliżone p-value (Stephens 1974)
    if A2_star < 0.2:
        p = 1 - math.exp(-13.436 + 101.14 * A2_star - 223.73 * A2_star**2)
    elif A2_star < 0.34:
        p = 1 - math.exp(-8.318 + 42.796 * A2_star - 59.938 * A2_star**2)
    elif A2_star < 0.6:
        p = math.exp(0.9177 - 4.279 * A2_star - 1.38 * A2_star**2)
    else:
        p = math.exp(1.2937 - 5.709 * A2_star + 0.0186 * A2_star**2)
    return A2_star, p


def mann_whitney_u(a, b):
    """Test Manna-Whitneya U (dwustronny) z poprawką na wiązania.
    H0: obie grupy pochodzą z tego samego rozkładu."""
    n1, n2 = len(a), len(b)
    all_vals = [(x, 1) for x in a] + [(x, 2) for x in b]
    all_vals.sort(key=lambda t: t[0])

    ranks = []
    i = 0
    while i < len(all_vals):
        j = i
        while j < len(all_vals) and all_vals[j][0] == all_vals[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2
        for k in range(i, j):
            ranks.append((all_vals[k][1], avg_rank))
        i = j

    r1 = sum(r for g, r in ranks if g == 1)
    u1 = r1 - n1 * (n1 + 1) / 2
    u2 = n1 * n2 - u1
    u = min(u1, u2)

    mean_u = n1 * n2 / 2
    N = n1 + n2
    counts = Counter(v for v, _ in all_vals)
    tie_corr = sum(c**3 - c for c in counts.values())
    std_u = math.sqrt((n1 * n2 / 12) * ((N + 1) - tie_corr / (N * (N - 1)))) if N > 1 else 1

    z = (u - mean_u) / std_u if std_u > 0 else 0
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return u, p


if __name__ == '__main__':
    data = load_csv('gc_content_results.csv')
    cats = ['INTERGENIC', 'INTRONIC', 'CDS']

    print('=' * 64)
    print('  1. TEST NORMALNOŚCI (Anderson-Darling) — GC na 5\'')
    print('=' * 64)
    print()
    print(f'{"Kategoria":<15} {"n":>6} {"A^2":>10} {"p-value":>12} {"Wniosek":<20}')
    print('-' * 64)
    for c in cats:
        A2, p = anderson_darling(data[c])
        if p is None:
            print(f'{c:<15} {len(data[c]):>6} {"n/d":>10} {"n/d":>12}')
            continue
        concl = 'NIE normalny' if p < 0.05 else 'normalny (n.s.)'
        print(f'{c:<15} {len(data[c]):>6} {A2:>10.2f} {p:>12.3e} {concl:<20}')
    print()
    print('Wniosek: jeśli p < 0.05, rozkład nie jest normalny -> test nieparametryczny.')
    print()

    print('=' * 64)
    print('  2. TEST MANN-WHITNEY U — GC na 5\'')
    print('=' * 64)
    print()
    for i in range(len(cats)):
        for j in range(i + 1, len(cats)):
            a, b = cats[i], cats[j]
            u, p = mann_whitney_u(data[a], data[b])
            sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'n.s.'))
            print(f'  {a:<12} vs {b:<12} n1={len(data[a])} n2={len(data[b])} U={u:.1f} p={p:.3e}  {sig}')
    print()
