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

## Wyniki analizy GC content (koniec 5')

| Kategoria | n | GC 5' (średnia) | GC całość (średnia) |
|-----------|-----|-----------------|---------------------|
| CDS       | 398 | **0.596**       | 0.524 |
| Intergenic| 941 | 0.518           | 0.476 |
| Intronic  | 1038| 0.509           | 0.468 |

**Test Shapiro-Wilk** (normalność): wszystkie kategorie p < 0.05 → rozkład nie normalny.

**Test Mann-Whitney U** (GC na 5'):

| Porównanie | p-value | Istotność |
|------------|---------|-----------|
| Intergenic vs Intronic | 0.084 | n.s. |
| Intergenic vs CDS      | <0.001 | *** |
| Intronic vs CDS        | <0.001 | *** |
