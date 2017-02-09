'''
genomediff.py

Part of parallyze
by Elizabeth Baird


Contains class and parsing function(s) for genomediff files.

'''

import sys
import os

import utils
from utils import complementary_base


class GenomeDiff(object):

    '''
    Base class for genomediffs. Attributes specific to mutation types
    are added after instantiation in the parsing function.
    '''

    def __init__(self, mut_id, mut_type, parent_ids, seq_id, position, fname):
        self.mut_type = mut_type
        self.mut_id = mut_id
        self.parent_ids = parent_ids
        self.seq_id = seq_id
        self.position = position
        self.fname = fname

    def __repr__(self):
        return '{0}: {1} {2}'.format(hash(self), self.mut_type, self.seq_id)


def parse_genomediff(gd_file, gb_record, genomediffs=None):
    '''
    Given a genomediff file, parse out the mutations and store them
    in GenomeDiff objects.

    If a dict for genomediffs is provided, add the mutations
    to that dict. Otherwise, create a new one.

    Returns a dictionary of all mutations keyed by mut_id.
    '''

    if genomediffs is None:
        genomediffs = {}
    else:
        assert isinstance(genomediffs, dict)

    with open(gd_file) as fp:
        for line in fp:
            # if line.startswith(i) for i in disregard_evidence:
            if line.startswith('#') or \
                    line.startswith('JC') or \
                    line.startswith('RA') or \
                    line.startswith('UN') or \
                    line.startswith('MC') or \
                    line.startswith('NOTE'):
                continue

            line = [token.strip() for token in line.split('\t')]

            '''
            The first columns are common to all GD's, so use them
            for the GenomeDiff object core
            '''

            mut_type = line[0]
            mut_id = line[1]
            parent_ids = line[2].split(',')
            seq_id = line[3]
            # NOTE: DOCUMENT THE ZERO-INDEXING OF POSITIONS IN GENOME!
            position = int(line[4]) - 1
            # which file (genome) did I come from?
            fname = gd_file

            gd = GenomeDiff(
                mut_id,
                mut_type,
                parent_ids,
                seq_id,
                position,
                fname)

            if mut_type == 'SNP':

                gd.new_base = line[5]
                assert gd.new_base in utils.dna
                '''
                Parse key/value pairs
                '''
                for pair in line[6:]:
                    key, _, value = pair.partition('=')
                    setattr(gd, key.strip(), value.strip())

                assert hasattr(gd, 'snp_type')
                assert hasattr(gd, 'gene_name')
                gd.gene_name = [gene.strip()
                                for gene in gd.gene_name.split('/')]
                assert hasattr(gd, 'locus_tag')
                gd.locus_tag = [tag.strip() for tag in gd.locus_tag.split('/')]
                assert hasattr(gd, 'gene_product')
                gd.gene_product = [product.strip()
                                   for product in gd.gene_product.split('/')]

                if gd.snp_type in ['nonsynonymous', 'synonymous']:

                    assert hasattr(gd, 'gene_strand')
                    if gd.gene_strand == '<':
                        gd.new_base = complementary_base(line[5])

                    assert hasattr(gd, 'codon_position')
                    assert hasattr(gd, 'codon_ref_seq')
                    gd.codon_position = int(gd.codon_position) - 1
                    gd.old_base = gd.codon_ref_seq[gd.codon_position]

                elif gd.snp_type in ['intergenic', 'pseudogene', 'noncoding']:

                    gd.old_base = gb_record.seq[gd.position]

                    if gd.new_base == gd.old_base:
                        print('WARNING: new base == old base' , \
                            gd.snp_type, ' problem with reference or indexing', \
                            file=sys.stderr)

                else:
                    print('ERROR: unrecognized SNP type ', snp_type, file=sys.stderr)

            elif mut_type == 'SUB':

                gd.size = int(line[5])
                gd.new_seq = line[6]

            elif mut_type == 'DEL':

                gd.size = int(line[5])

            elif mut_type == 'INS':

                gd.new_seq = line[5]

            elif mut_type == 'MOB':

                gd.repeat_name = line[5]
                gd.strand = line[6]
                gd.duplication_size = int(line[7])

            elif mut_type == 'AMP':

                gd.size = int(line[5])
                gd.new_copy_number = int(line[6])

            elif mut_type == 'CON':

                gd.size = int(line[5])
                gd.region = line[6]

            elif mut_type == 'INV':
                gd.size = int(line[5])

            # Get a unique key for addressing this mutation  ## NOTE: I don't
            # like this hashing but leave for now...what is a better key?
            key = hash(gd)
            # Store the line for this mutation
            gd.line = os.path.basename(gd_file)

            genomediffs[key] = gd

    return genomediffs
