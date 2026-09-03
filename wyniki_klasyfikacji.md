# Wyniki klasyfikacji retrogenów

## Tabela zbiorcza

| Lista | Kategoria | Liczba | ENCODE | % | Kryterium przyjęcia |
|-------|-----------|--------|--------|---|---------------------|
| **1** | Międzygenowe | **5 205** | 941 | 18.1% | 0% overlapu z transkryptem + self-hity |
| 1a | — oryginalne | 1 949 | 339 | 14.2% | zero overlapu z jakimkolwiek transkryptem ncbiRefSeq |
| 1b | — self-hit | 3 256 | 820 | 25.2% | w transkrypcie, ale tylko własnym (brak gospodarza) |
| **2** | Intronowe | **3 191** | 1 055 | 33.1% | ≥50% w transkrypcie + 0% overlapu z eksonem |
| **3** | CDS 1 egzon | **768** | 398 | 51.8% | ≥50% w CDS transkryptu z dokładnie 1 egzonem kodującym |
| **4** | Reszta | **5 710** | 1 669 | 29.2% | dotyka eksonu (jakikolwiek overlap) — stany graniczne |
| | **Suma** | **14 874** | | | |

## Uwagi

- **Self-hit (1b)** — podzbiór międzygenowych: retrogeny które w adnotacji ncbiRefSeq
  znalazły tylko same siebie (własny pseudogen), a nie „siedzą" w innym genie.
  Wyróżnione w osobnym `list5_selfhit.bed`.
- **Priorytet klasyfikacji**: CDS 1 egzon → intron → reszta → intergenic.
- **% ENCODE** — odsetek retrogenów danej kategorii z potwierdzeniem w długich odczytach
  ENCODE (po liftoverze hg38→hs1), próg `-f 0.1 -F 0.1`.

## Wyniki analizy GC content (podział na połowy, jak Mordstein 2020)

| Kategoria | n | GC 1. połowa (5') | GC 2. połowa (3') | Różnica |
|-----------|-----|-------------------|-------------------|---------|
| CDS       | 398 | **57.0%**         | 47.7%             | +9.3 p.p. |
| Intergenic| 941 | 49.3%             | 45.9%             | +3.4 p.p. |
| Intronic  | 1038| 48.3%             | 45.3%             | +3.0 p.p. |

**Test Wilcoxona dla par** (pierwsza vs druga połowa): wszystkie kategorie
istotnie wyższe GC na 5' (p < 0.001).

**Test Shapiro-Wilk** (normalność GC 1. połowy): wszystkie kategorie p < 0.05
→ rozkład nie normalny → test nieparametryczny.

**Test Mann-Whitney U** (GC 1. połowy między kategoriami, z korektą Bonferroniego):

| Porównanie | p-value | Istotność |
|------------|---------|-----------|
| Intergenic vs Intronic | 0.014 | * |
| Intergenic vs CDS      | 1.4e-41 | *** |
| Intronic vs CDS        | 8.9e-55 | *** |
