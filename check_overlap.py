import pandas as pd

#czwarte kolumny (nazwy retrokopii) z każdego pliku jako zestawy (set)
def get_names(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, usecols=[3])
    return set(df[3].tolist())

cds_retrogenes = get_names('retrogenes_in_single_exon_cds.bed')
intergenic = get_names('intergenic_retrogenes.bed')
intronic = get_names('intronic_retrogenes_genes.bed')

# szukamy części wspólnych (iloczyn zbiorów, powinny być zbiorem rozłącznym)
overlap_cds_intergenic = cds_retrogenes.intersection(intergenic)
overlap_cds_intronic = cds_retrogenes.intersection(intronic)
overlap_intergenic_intronic = intergenic.intersection(intronic)

#Wyświetlamy wyniki
print("--- WYNIKI KONTROLI ROZŁĄCZNOŚCI ---")

if not overlap_cds_intergenic and not overlap_cds_intronic and not overlap_intergenic_intronic:
    print("Sukces! Pliki są całkowicie rozłączne. Żaden retrogen się nie powtarza.")
else:
    if overlap_cds_intergenic:
        print(f"Duble między Single CDS a Intergenic ({len(overlap_cds_intergenic)}):")
        print(overlap_cds_intergenic)
        
    if overlap_cds_intronic:
        print(f"Duble między Single CDS a Intronic ({len(overlap_cds_intronic)}):")
        print(overlap_cds_intronic)
        
    if overlap_intergenic_intronic:
        print(f"Duble między Intergenic a Intronic ({len(overlap_intergenic_intronic)}):")
        print(overlap_intergenic_intronic)
