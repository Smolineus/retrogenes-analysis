#!/usr/bin/env python3
"""
Porównanie statystyczne GC content (5' koniec) między kategoriami.
Test Mann-Whitney U (nieparametryczny, dane nie mają rozkładu normalnego).
"""

import csv
from collections import defaultdict

def load_csv(path):
    data = defaultdict(list)
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            data[row['kategoria']].append(float(row['gc_5prime']))
    return data

def mann_whitney_u(a, b):
    """Test Manna-Whitneya U (dwustronny), zwraca statystykę U i p-value (przybliżenie normalne)."""
    import math
    n1, n2 = len(a), len(b)
    # rank wszystkich wartości
    all_vals = [(x, 1) for x in a] + [(x, 2) for x in b]
    all_vals.sort(key=lambda t: t[0])

    # rankingi z poprawką na wiązania
    ranks = []
    i = 0
    while i < len(all_vals):
        j = i
        while j < len(all_vals) and all_vals[j][0] == all_vals[i][0]:
            j += 1
        avg_rank = (i + j + 1) / 2
        for k in range(i, j):
            ranks.append((all_vals[k][0], all_vals[k][1], avg_rank))
        i = j

    r1 = sum(r for _, g, r in ranks if g == 1)
    u1 = r1 - n1 * (n1 + 1) / 2
    u2 = n1 * n2 - u1
    u = min(u1, u2)

    # przybliżenie normalne
    mean_u = n1 * n2 / 2
    # poprawka na wiązania
    # liczba wiązań
    from collections import Counter
    counts = Counter(v for v, g in all_vals)
    tie_corr = 0
    N = n1 + n2
    for c in counts.values():
        tie_corr += c**3 - c
    std_u = math.sqrt((n1 * n2 / 12) * ((N + 1) - tie_corr / (N * (N - 1)))) if N > 1 else 1

    z = (u - mean_u) / std_u if std_u > 0 else 0
    # p-value dwustronne (przybliżenie)
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return u, p

if __name__ == '__main__':
    data = load_csv('gc_content_results.csv')
    cats = ['INTERGENIC', 'INTRONIC', 'CDS']
    print('=' * 60)
    print('  TEST MANN-WHITNEY U — GC content na końcu 5\'')
    print('=' * 60)
    print()
    for i in range(len(cats)):
        for j in range(i+1, len(cats)):
            a, b = cats[i], cats[j]
            u, p = mann_whitney_u(data[a], data[b])
            sig = '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else 'n.s.'))
            print(f'  {a} vs {b}:')
            print(f'    n1={len(data[a])}, n2={len(data[b])}, U={u:.1f}, p={p:.3e}  {sig}')
    print()
