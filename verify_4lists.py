#!/usr/bin/env python3
import sys

FILES = {
    'list1_intergenic.bed':       'Intergenowe',
    'list2_intronic.bed':         'Intronowe',
    'list3_cds_one_exon.bed':     'CDS 1 egzon',
    'list4_exonic_multigene.bed': 'Eksonowe multi-gen',
}

def get_names(filepath):
    names = set()
    try:
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if line:
                    names.add(line.split('\t')[3])
    except FileNotFoundError:
        print(f"  NIE ZNALEZIONO: {filepath}")
    return names

def main():
    all_sets = {}
    for fpath, label in FILES.items():
        all_sets[label] = get_names(fpath)

    original = get_names('human_retrocopies_merged.bed')

    print("=" * 55)
    print("  LICZNOŚCI")
    print("=" * 55)
    for label, names in all_sets.items():
        print(f"  {label:25s} {len(names):>6d}")
    all_assigned = set()
    for names in all_sets.values():
        all_assigned |= names
    print(f"  {'-' * 35}")
    print(f"  Suma (4 listy):          {len(all_assigned):>6d}")
    print(f"  Oryginalny plik:         {len(original):>6d}")

    print()
    print("=" * 55)
    print("  KONTROLA ROZŁĄCZNOŚCI")
    print("=" * 55)
    labels = list(all_sets.keys())
    disjoint = True
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            overlap = all_sets[labels[i]] & all_sets[labels[j]]
            if overlap:
                disjoint = False
                print(f"  KONFLIKT: {labels[i]} & {labels[j]}: {len(overlap)}")
            else:
                print(f"  OK:       {labels[i]} & {labels[j]}: rozłączne")

    print()
    print("=" * 55)
    print("  BRAKUJĄCE / NADMIAROWE")
    print("=" * 55)
    missing = original - all_assigned
    extra = all_assigned - original

    if missing:
        print(f"  Brakuje:  {len(missing)} retrogenów")
        for name in sorted(list(missing)[:10]):
            print(f"    - {name}")
        if len(missing) > 10:
            print(f"    ... i {len(missing) - 10} więcej")
    else:
        print(f"  Brakuje:  0")

    if extra:
        print(f"  Nadmiar:  {len(extra)} retrogenów (spoza oryginalnego pliku!)")
        for name in sorted(list(extra)[:10]):
            print(f"    + {name}")
    else:
        print(f"  Nadmiar:  0")

    print()
    if not missing and not extra and disjoint:
        print("  SUKCES: Wszystkie retrogeny sklasyfikowane, listy rozłączne.")
    elif not missing and disjoint:
        print("  UWAGA: Listy pokrywają wszystkie retrogeny, są rozłączne, ale mają nadmiarowe wpisy.")
    elif not disjoint:
        print("  UWAGA: Listy NIE są rozłączne!")
    elif missing:
        print(f"  UWAGA: Brakuje {len(missing)} retrogenów!")

if __name__ == '__main__':
    sys.exit(main() or 0)
