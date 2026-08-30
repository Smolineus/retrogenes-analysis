import pandas as pd
import gzip

# nazywanie kolumn w gtf na podstawie norm
def get_one_CDS_transcript(gtf_path, output_bed):
    cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
   
    print('Wczytywanie GTF...')
    # wczytywanie lini z CDS 
    all_CDS = [] 
    with gzip.open(gtf_path, 'rt') as f: # rt - read text -> tryb odczytu tekstowego
        for line in f: 
            fields = line.strip().split('\t')
            if fields[2] in ['CDS']:
                all_CDS.append(fields)

    df = pd.DataFrame(all_CDS, columns=cols)
    
    # Wyciąganie gene_id oraz transcript_id
    def get_gene_and_transcript(attr):
        # unknown na wszelki wypadek jakby gdzieś nie było podane, żeby błędu nie wyrzuciło
        g_id = "unknown_gene"
        t_id = "unknown_transcript"
        
        for item in attr.split(';'): # bo ostatnia kolumna oddzielona ;
            item = item.strip() # usuwanie spacji z brzegów
            
            if item.startswith('gene_id'):
                g_id = item.split('"')[1]
                
            elif item.startswith('transcript_id'):
                t_id = item.split('"')[1]
                
        # Zwracamy połączony napis w formacie: gene_id: X; transcript_id: Y
        return f'gene_id: {g_id}; transcript_id: {t_id}'

    # Funkcja pomocnicza do wyciągania SAMEGO transcript_id (potrzebna do grupowania)
    def get_just_transcript_id(attr):
        for item in attr.split(';'):
            if 'transcript_id' in item:
                return item.split('"')[1]
        return "unknown_transcript"

    # Tworzymy nowe kolumny
    df['full_id'] = df['attribute'].apply(get_gene_and_transcript)
    df['transcript_id'] = df['attribute'].apply(get_just_transcript_id)

    # liczenie ile dany transcript id ma przypisanych CDS
    cds_counts = df.groupby('transcript_id').size()
    single_cds_id = cds_counts[cds_counts == 1].index.tolist()

    print(f'Znaleziono {len(single_cds_id)}')

    # filtrowanie
    final_cds_df = df[df['transcript_id'].isin(single_cds_id)].copy()

    # konwersja na format BED (0-based start)
    result_bed = final_cds_df[['seqname', 'start', 'end', 'full_id', 'strand']].copy()
    result_bed['start'] = result_bed['start'].astype(int) - 1 
    result_bed['end'] = result_bed['end'].astype(int)

    # zapis do pliku
    result_bed.to_csv(output_bed, sep='\t', header=False, index=False)
    print(f'Zapisano do pliku: {output_bed}')


# Kod główny (wywołanie funkcji) - od samego brzegu pliku
input_gtf = 'hs1.ncbiRefSeq.bigZip.gtf.gz'
output_name = 'single_cds.bed'

get_one_CDS_transcript(input_gtf, output_name)
