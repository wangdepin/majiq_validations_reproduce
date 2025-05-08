from majiq.src.gff import parse_gff3
import argparse
import sys

accepted_transcripts = ['mRNA', 'transcript', 'lnc_RNA', 'miRNA', 'ncRNA',
                                  'rRNA', 'scRNA', 'snRNA', 'snoRNA', 'tRNA', 'pseudogenic_transcript',
                                  'C_gene_segment', 'D_gene_segment', 'J_gene_segment',
                                  'V_gene_segment', 'unconfirmed_transcript', 'three_prime_overlapping_ncrna']
transcript_id_keys = 'ID'
accepted_genes = ['gene', 'ncRNA_gene', 'pseudogene', 'ncRNA_gene', 'bidirectional_promoter_lncRNA']
gene_name_keys = ['Name', 'gene_name']
gene_id_keys = ['ID', 'gene_id']


def read_gff(filename):

    all_genes = {}
    transcript_len = {}
    for record in parse_gff3(filename):

        if record.strand is None or record.seqid is None:
            continue

        # print ("### ", record)
        start = record.start
        end = record.end
        if record.type in accepted_genes:
            for gname_k in gene_name_keys:
                try:
                    gene_name = record.attributes[gname_k]
                    break
                except KeyError:
                    continue
            else:
                print("Error, Gene doesn't contain one of the Name attribute  information values: "
                      "%s" % gene_name_keys)
            for gid_k in gene_id_keys:
                try:
                    gene_id = record.attributes[gid_k]
                    break
                except KeyError:
                    continue
            else:
                print("Error, Gene doesn't contain one of the ID attribute information values: "
                      "%s" % gene_id_keys)
            if gene_id in all_genes:
                raise RuntimeError('Two Different Genes with the same name %s' % gene_name)
            all_genes[gene_id] = {}

        elif record.type in accepted_transcripts:
            if transcript_id_keys not in record.attributes or 'Parent' not in record.attributes:
                print("Error, Transcript doesn't contain one of the ID or parent attributes"
                      "information values: %s" % transcript_id_keys)
                continue
            transcript_name = record.attributes[transcript_id_keys]
            parent = record.attributes['Parent']
            if gene_id not in all_genes:
                print("Error, incorrect gff. mRNA %s doesn't have valid gene %s" % (transcript_name, parent))
                continue
            transcript_len[record.attributes['ID']] = 0

        elif record.type == 'exon':
            parent_tx_id = record.attributes['Parent']
            try:
                transcript_len[parent_tx_id] +=  (end-start)

            except KeyError:
                print("Error, incorrect gff. exon "
                      "doesn't have valid mRNA %s" % parent_tx_id)

    return transcript_len


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('ttable', action='store')
    parser.add_argument('db', action='store')
    parser.add_argument('output_file', action='store')
    args = parser.parse_args()

    tnscpt_len = read_gff(args.db)

    with open(args.ttable) as fp , open(args.output_file, 'w+') as ofp:
        ofp.write('#This file has been generated with the following command:\n')
        ofp.write('# %s\n' % ' '.join(sys.argv))
        ofp.write('#SIM_ID\t#TXT_ID\tGENE_ID\tTXT_LENGTH\n')
        for xx in fp.readlines():
            if xx.startswith('#'): continue
            tab = xx.strip().split()
            try:
                ln = tnscpt_len[tab[1]]
            except KeyError:
                print ('ERROR TXT', tab[1])
                continue
            ofp.write('%s\t%s\n' % ('\t'.join(tab), ln))
