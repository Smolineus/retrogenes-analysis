#!/usr/bin/env python3
"""
Generowanie wykresów do prezentacji — wyniki klasyfikacji, walidacji ENCODE
i analizy GC. Wykresy zapisywane jako PNG w katalogu charts/.
"""

import os
import csv
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')  # backend bez okna (do zapisu plików)
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

os.makedirs('charts', exist_ok=True)

# Styl
plt.rcParams.update({
    'figure.dpi': 150,
    'font.size': 11,
    'axes.spines.top': False,
    'axes.spines.right': False,
})

# Paleta kolorów (spójna we wszystkich wykresach)
# Międzygenowe -> błękitny, Intronowe -> różowy, CDS -> fioletowy
BLUE = '#4FC3F7'
PINK = '#F06292'
PURPLE = '#BA68C8'

CATS = ['INTERGENIC', 'INTRONIC', 'CDS']
LABELS = ['Międzygenowe', 'Intronowe', 'CDS 1 egzon']
COLORS = [BLUE, PINK, PURPLE]


# ============ funkcje pomocnicze ============

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


def load_halves():
    data = defaultdict(lambda: {'gc_first': [], 'gc_second': [], 'lengths': []})
    with open('fasta/gc_halves_results.csv') as f:
        reader = csv.DictReader(f)
        for row in reader:
            k = row['kategoria']
            data[k]['gc_first'].append(float(row['gc_pierwsza_polowa']) * 100)
            data[k]['gc_second'].append(float(row['gc_druga_polowa']) * 100)
            data[k]['lengths'].append(int(row['dlugosc']))
    return data


def normalized_gc_profile(seq, n_bins=20):
    """GC content w n_bins znormalizowanych przedziałach (0% = 5', 100% = 3')."""
    if len(seq) < n_bins:
        return None
    profile = []
    for i in range(n_bins):
        start = int(i * len(seq) / n_bins)
        end = int((i + 1) * len(seq) / n_bins)
        chunk = seq[start:end]
        gc = sum(1 for c in chunk if c in 'GC') / len(chunk)
        profile.append(gc)
    return profile


data = load_halves()

# ============ WYKRES 1: Klasyfikacja retrogenów (3 czyste kategorie) ============
cats_clean = ['CDS 1 egzon', 'Intronowe', 'Międzygenowe']
counts_clean = [768, 3191, 5205]
colors_clean = [PURPLE, PINK, BLUE]

fig, ax = plt.subplots(figsize=(8, 4.5))
bars = ax.barh(cats_clean, counts_clean, color=colors_clean, height=0.6)
for bar, c in zip(bars, counts_clean):
    ax.text(bar.get_width() + 150, bar.get_y() + bar.get_height()/2,
            f'{c:,}'.replace(',', ' '), va='center', fontweight='bold')
ax.set_xlabel('Liczba retrogenów')
ax.set_title('Klasyfikacja retrogenów — kategorie czyste')
ax.set_xlim(0, 6000)
ax.invert_yaxis()  # największa na górze
plt.tight_layout()
plt.savefig('charts/klasyfikacja.png', bbox_inches='tight')
plt.close()

# ============ WYKRES 2: Wsparcie ENCODE (%) ============
categories = LABELS
enc_pct = [18.1, 33.1, 51.8]

fig, ax = plt.subplots(figsize=(8, 5))
bars = ax.bar(categories, enc_pct, color=COLORS)
for bar, p in zip(bars, enc_pct):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
            f'{p:.1f}%', ha='center', va='bottom', fontweight='bold')
ax.set_ylabel('Odsetek retrogenów z wsparciem ENCODE [%]')
ax.set_title('Walidacja ekspresji — wsparcie długimi odczytami ENCODE')
ax.set_ylim(0, 60)
plt.tight_layout()
plt.savefig('charts/wsparcie_encode.png', bbox_inches='tight')
plt.close()

# ============ WYKRES 3: GC — połówki (pierwsza vs druga) ============
gc_first_means = [np.mean(data[c]['gc_first']) for c in CATS]
gc_second_means = [np.mean(data[c]['gc_second']) for c in CATS]

x = np.arange(len(LABELS))
width = 0.35

fig, ax = plt.subplots(figsize=(8, 5))
b1 = ax.bar(x - width/2, gc_first_means, width, label='Pierwsza połowa (5\')', color=BLUE)
b2 = ax.bar(x + width/2, gc_second_means, width, label='Druga połowa (3\')', color=PINK)
for b in b1:
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.5,
            f'{b.get_height():.1f}%', ha='center', va='bottom', fontsize=9)
for b in b2:
    ax.text(b.get_x() + b.get_width()/2, b.get_height() + 0.5,
            f'{b.get_height():.1f}%', ha='center', va='bottom', fontsize=9)
ax.set_ylabel('Zawartość GC [%]')
ax.set_title('Zawartość GC — pierwsza vs druga połowa transkryptu')
ax.set_xticks(x)
ax.set_xticklabels(LABELS)
ax.legend()
ax.set_ylim(0, 65)
plt.tight_layout()
plt.savefig('charts/gc_polowy.png', bbox_inches='tight')
plt.close()

# ============ WYKRES 4: Boxplot GC 5' między kategoriami ============
gc_first_by_cat = [data[c]['gc_first'] for c in CATS]

fig, ax = plt.subplots(figsize=(8, 5))
bp = ax.boxplot(gc_first_by_cat, tick_labels=LABELS, patch_artist=True,
                widths=0.5, showfliers=False, medianprops=dict(color='black', linewidth=2))
for patch, color in zip(bp['boxes'], COLORS):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

for i, vals in enumerate(gc_first_by_cat):
    med = np.median(vals)
    ax.text(i + 1, med + 0.5, f'med={med:.1f}%', ha='center', va='bottom', fontsize=9)

ax.set_ylabel('Zawartość GC (pierwsza połowa) [%]')
ax.set_title('Rozkład GC na 5\' — porównanie kategorii')
ax.set_ylim(20, 95)

plt.tight_layout()
plt.savefig('charts/boxplot_gc5.png', bbox_inches='tight')
plt.close()

# ============ WYKRES 5: Boxplot długości sekwencji ============
lens_by_cat = [data[c]['lengths'] for c in CATS]

fig, ax = plt.subplots(figsize=(8, 5))
bp = ax.boxplot(lens_by_cat, tick_labels=LABELS, patch_artist=True,
                widths=0.5, showfliers=False, medianprops=dict(color='black', linewidth=2))
for patch, color in zip(bp['boxes'], COLORS):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)

for i, vals in enumerate(lens_by_cat):
    med = np.median(vals)
    ax.text(i + 1, med, f'med={med:.0f} nt', ha='center', va='bottom', fontsize=9)

ax.set_yscale('log')
ax.set_ylabel('Długość transkryptu [nt] (skala log)')
ax.set_title('Długość transkryptów — porównanie kategorii')
plt.tight_layout()
plt.savefig('charts/histogram_dlugosci.png', bbox_inches='tight')
plt.close()

# ============ WYKRES 6: Gradient GC 5'->3' (znormalizowany) ============
N_BINS = 20  # 5% długości na bin

fasta_files = {
    'INTERGENIC': 'fasta/intergenic_transcripts.fa',
    'INTRONIC': 'fasta/intronic_transcripts.fa',
    'CDS': 'fasta/cds_transcripts.fa',
}

fig, ax = plt.subplots(figsize=(8, 5))
x_axis = np.linspace(0, 100, N_BINS)

for c, color, label in zip(CATS, COLORS, LABELS):
    seqs = read_fasta(fasta_files[c])
    profiles = []
    for seq in seqs.values():
        p = normalized_gc_profile(seq)
        if p is not None:
            profiles.append(p)
    mean_profile = np.mean(profiles, axis=0) * 100
    ax.plot(x_axis, mean_profile, color=color, label=label, lw=2.5)

ax.set_xlabel('Pozycja w transkrypcie [%] (0% = 5\', 100% = 3\')')
ax.set_ylabel('Zawartość GC [%]')
ax.set_title('Gradient GC wzdłuż transkryptu (znormalizowany)')
ax.set_ylim(42, 62)
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig('charts/gradient_gc.png', bbox_inches='tight')
plt.close()

print('Wykresy zapisane w charts/:')
for f in ['klasyfikacja.png', 'wsparcie_encode.png', 'gc_polowy.png',
          'boxplot_gc5.png', 'histogram_dlugosci.png', 'gradient_gc.png']:
    print(f'  {f}')
