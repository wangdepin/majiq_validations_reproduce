# This script reads a gff3 file and extracts all the annotated introns including gene id transcript Id and intron coordinates
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

gene_coords = {}
gene2chrom = {}

def check_overlap(gid, start, end, glist):
    r = False
    for vv in glist:
        if vv[0] < end and vv[1]> start and vv[2] != gid :
            return True
    return False

def read_gff(filename):

    trcpt_id_dict = {}
    exon_dict = {}
    all_genes = {}
    junc_set = set()
    for record in parse_gff3(filename):

        if record.strand is None or record.seqid is None:
            continue

        # print ("### ", record)
        chrom = record.seqid
        strand = record.strand
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

            exon_dict[gene_id] = []
            all_genes[gene_id] = {}

            k = '%s:%s' %(chrom, strand)
            gene2chrom[gene_id] = k
            try:
                gene_coords[k].append((start, end, gene_id))
            except KeyError:
                gene_coords[k] = [(start, end, gene_id)]

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

            trcpt_id_dict[record.attributes['ID']] = [parent, []]

        elif record.type == 'exon':
            parent_tx_id = record.attributes['Parent']
            try:
                gene_id = trcpt_id_dict[parent_tx_id][0]
                exon_dict[gene_id].append((start, True,  parent_tx_id))
                exon_dict[gene_id].append((end,   False, parent_tx_id))
                trcpt_id_dict[parent_tx_id][1].append((start, end))

            except KeyError:
                print("Error, incorrect gff. exon "
                      "doesn't have valid mRNA %s" % parent_tx_id)

    for parent_tx_id, (gene_id, coord_list) in trcpt_id_dict.items():
        last_ss = -1
        coord_list.sort(key=lambda x: (x[0], x[1]))
        if len(coord_list) == 0: continue
        for xx, yy in coord_list:
            key = '%s-%s' % (last_ss, xx)
            junc_set.add(key)
            last_ss = yy


    merge_exons(exon_dict, all_genes, junc_set)
    return all_genes


def merge_exons(exon_dict, all_genes, junc_set):

    ltxt = set()
    for gne_id, ex_list in exon_dict.items():
        ex_list.sort(key=lambda x:(x[0], -x[1]))
        ex_start = -1
        ex_end = -1
        nopen = 0

        for coord, is_start, tid in ex_list:
            # print( gne_id )
            # if gne_id == 'ENSG00000254413' :
            #     print(gne_id, coord, is_start, tid, nopen, ltxt)

            if is_start:
                if ex_end != -1:
                    if nopen > 0 and (ex_end+4) < (coord-1):
                        for transcript_id in ltxt:
                            ir_key = "%s-%s" % (ex_end+1, coord-1)
                            j_c = "%s-%s" % (ex_end, coord)
                            if j_c not in junc_set:
                                continue

                            if not check_overlap(gne_id, ex_end+1, coord-1, gene_coords[gene2chrom[gne_id]]):
                                try:
                                    all_genes[gne_id][ir_key].append(transcript_id)
                                except KeyError:
                                    all_genes[gne_id][ir_key] = [transcript_id]

#                    else:
                        # all_genes[gne_id].create_annot_intron(ex_end+1, coord-1)
                        #tlist.append([ex_end+1, coord-1, 1, IR_TYPE])
                    ex_end = -1
                    nopen = 0
                    ex_start = coord

                ex_start = coord if ex_start == -1 or coord < ex_start else ex_start
                nopen += 1
                # print("ADD", tid)
                ltxt.add(tid)

            else:
                nopen -= 1
                ex_end = coord if coord > ex_end else ex_end
                # print("REMOVE", tid)
                ltxt.remove(tid)


def extrac_introns(db_file, args):
    all_gene = read_gff(db_file)
    if (args.inclusion_file):
        pass

    with open(args.output_file, 'w+') as ofp:
        ofp.write('#This file has been generated with the following command:\n')
        ofp.write('# %s\n' % ' '.join(sys.argv))
        ofp.write('#GENE_ID\t#IR_COORDINATES\tTRANSCRIPT_IDS\n' )
        for gg, rr in all_gene.items():
            for kk, lv in rr.items():
                ofp.write("%s\t%s\t%s\n" % (gg, kk, ";".join(lv)))





if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('db', action='store')
    parser.add_argument('output_file', action='store')
    parser.add_argument('-i', '--inclusion_file', action='store')
    args = parser.parse_args()

    extrac_introns(args.db, args)
