'''
gene.py

Defines a Gene object, given a feature object from Biopython's
SeqIO module.
'''


class Gene(object):

    def __init__(self, feature):
        assert feature.type == 'CDS'
        self.locus_tag = feature.qualifiers['locus_tag'][0]
        try:
            self.name = feature.qualifiers['gene'][0]
        except KeyError:
            self.name = None
        self.start = int(feature.location.start)
        self.end = int(feature.location.end)
        self.A = 0
        self.G = 0
        self.C = 0
        self.T = 0

    def count_bases(self, refseq):
        '''
        Given an iterator defining a genomic sequence,
        count the number of each nucleotide type within the
        coordinates of this gene.
        '''
        for base in refseq[self.start:self.end]:
            if base in ['A', 'a']:
                self.A += 1
            elif base in ['T', 't']:
                self.T += 1
            elif base in ['C', 'c']:
                self.C += 1
            elif base in ['G', 'g']:
                self.G += 1
            else:
                pass

    def __str__(self):
        return '{} ({}): ({},{}) A={} G={} C={} T={}'.format(
            self.name, self.locus_tag, self.start, self.end,
            self.A, self.G, self.C, self.T)


def get_genecoordinates(record):
    geneinfo = {}
    total_bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'total': 0}
    refseq = list(record.seq)
    for feature in record.features:
        if feature.type == 'CDS':
            this_gene = Gene(feature)
            this_gene.count_bases(refseq)
            geneinfo[this_gene.locus_tag] = this_gene
            total_bases['A'] += this_gene.A
            total_bases['G'] += this_gene.G
            total_bases['C'] += this_gene.C
            total_bases['T'] += this_gene.T
    total_bases['total'] = total_bases['A'] + \
        total_bases['G'] + total_bases['C'] + total_bases['T']
    print('\n', 'total A, G, C, T, TOTAL', total_bases['A'], \
        total_bases['G'], total_bases['C'], total_bases['T'], \
        total_bases['total'], sep=" ")
    # print('\n', 'Reference genome list (1st 10): gene name, locus_tag, \
        # start, stop, A, G, C, T', '\n', geneinfo[:10], sep=' ')
    # geneinfo is now unhashable type: b/c of data restructure?
    return geneinfo, total_bases
