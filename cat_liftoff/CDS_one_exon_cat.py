#!/usr/bin/env python3
import gzip

def get_one_CDS_transcript(gtf_path, output_bed):
    cds_by_transcript = {}
    cds_lines = []
    transcript_chr = {}

    print('Wczytywanie GTF...')
    with gzip.open(gtf_path, 'rt') as f:
        for line in f:
            fields = line.strip().split('\t')
            if fields[2] == 'CDS':
                attr = fields[8]
                tid = 'unknown'
                for item in attr.split(';'):
                    item = item.strip()
                    if item.startswith('transcript_id'):
                        tid = item.split('"')[1]
                        break
                cds_by_transcript[tid] = cds_by_transcript.get(tid, 0) + 1
                cds_lines.append(fields)
                transcript_chr[tid] = fields[0]

    single_cds_ids = {tid for tid, count in cds_by_transcript.items() if count == 1}
    print(f'Znaleziono {len(single_cds_ids)} transkryptow z 1 CDS')

    seen_keys = set()
    with open(output_bed, 'w') as out:
        for fields in cds_lines:
            attr = fields[8]
            tid = 'unknown'
            gid = 'unknown'
            for item in attr.split(';'):
                item = item.strip()
                if item.startswith('transcript_id'):
                    tid = item.split('"')[1]
                elif item.startswith('gene_id'):
                    gid = item.split('"')[1]

            if tid not in single_cds_ids:
                continue

            full_id = f'gene_id: {gid}; transcript_id: {tid}'
            seqname = fields[0]
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]

            key = (seqname, start, end, full_id, strand)
            if key in seen_keys:
                continue
            seen_keys.add(key)

            out.write(f'{seqname}\t{start}\t{end}\t{full_id}\t{strand}\n')

    print(f'Zapisano do pliku: {output_bed}')


input_gtf = 'catLiftOffGenesV1.gtf.gz'
output_name = 'single_cds_cat.bed'
get_one_CDS_transcript(input_gtf, output_name)
