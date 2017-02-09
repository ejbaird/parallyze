# parallyze_routines.py

import numpy as np
import matplotlib.pyplot as plt
import operator
import time
import os
import shutil

from genomediff import parse_genomediff, GenomeDiff
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.Applications import RaxmlCommandline
from Bio import AlignIO
from Bio import Phylo

import random
import utils
import itertools
from pprint import pprint
from functools import reduce


def fileName(args):
    timestr = time.strftime('%Y_%m_%d_%H_%M_%S')
    if args.fname is None:
        out_fn = timestr+'.tsv'
    else:
        out_fn = args.fname+timestr+'.tsv'

    return out_fn


def calculate_dNdS(genomediffs):
    '''input: the output from parsed gd files
    goal: count the SNP mutation types for dN/dS, intergenic, etc.
    output: dict of dict of a count'''

    counts = {}
    dNtotal = 0
    dStotal = 0
    dN1 = 0
    dS1 = 0
    dN2 = 0
    dS2 = 0
    dN3plus = 0
    dS3plus = 0
    loci1 = 0
    loci2 = 0
    loci3plus = 0

    for mut_id, gd in genomediffs.items():
        if gd.mut_type == 'SNP' and \
                gd.snp_type in ['synonymous', 'nonsynonymous']:
                # These all only have one locus_tag, and we can't
                # use a list as a key, so just get the value
            locus_tag = gd.locus_tag[0]

            if locus_tag in counts:
                dN, dS = counts[locus_tag]
            else:
                dN, dS = (0, 0)

            if gd.snp_type == 'nonsynonymous':
                dN += 1
                dNtotal += 1
            else:
                dS += 1
                dStotal += 1

            counts[locus_tag] = (dN, dS)

    for locusTag, dnds in counts.items():
        # # # if sum of tuple of dN and dS per locus tag =1,...etc.
        if dnds[0] + dnds[1] == 1:
            dN1 += dnds[0]
            dS1 += dnds[1]
            loci1 += 1
        elif dnds[0] + dnds[1] == 2:
            dN2 += dnds[0]
            dS2 += dnds[1]
            loci2 += 1
        elif dnds[0] + dnds[1] >= 3:
            dN3plus += dnds[0]
            dS3plus += dnds[1]
            loci3plus += 1

    try:
        dNdS1 = float(dN1) / float(dS1)
    except ZeroDivisionError:
        dNdS1 = dN1
        print("dS for singly mutating loci is 0. dN=", dN1)
        # pseudocount: assume dS=1 to allow for math

    try:
        dNdS2 = float(dN2) / float(dS2)
    except ZeroDivisionError:
        dNdS2 = dN2
        print("dS for doubly mutating loci is 0. dN=", dN2)

    try:
        dNdS3plus = float(dN3plus) / float(dS3plus)
    except ZeroDivisionError:
        print("dS for triply or more mutating loci is 0. dN=", dN3plus)
        dNdS3plus = dN3plus

    print("loci with 1 mutation:", loci1)
    print("loci with 2 mutations:", loci2)
    print("loci with 3 or more mutations:", loci3plus)

    return counts, dNtotal, dStotal, dNdS1, dNdS2, dNdS3plus


def dndsPlot(dndsGraphingDict):
    pass


# NOTE: Updated for refactor
def snpcount(genomediffs, lines, snp_types):
    '''input: 'mutations' i.e., parsed gd files
    returns dictionary of gdfiles, each containing a matrix
    of SNP mutations - to and from base'''

    file_matrices = {}
    snp_type_counts = {snp_types}
    print(snp_type_counts)
    for line in lines:
        file_matrices[line] = np.zeros((4, 4), dtype=int)

    for key, gd in genomediffs.items():
        if gd.mut_type == 'SNP' and gd.snp_type in snp_types:
            old_base = utils.seq_to_int(gd.old_base)
            new_base = utils.seq_to_int(gd.new_base)
            line = gd.line
            file_matrices[line][old_base, new_base] += 1

            # do i need to set each snp type count to 0?
            snp_type_counts[gd.snp_type] += 1

    for line, mat in file_matrices.items():
        print('file:', line)
        print('from (row) / to (column) :', '\n', mat, '\n')

    return file_matrices, snp_type_counts  # updated to old procs?


# NOTE: Updated for refactor
def mutated_lines_per_gene(genomediffs, snp_types):
    '''input: mutations, from parse_gdfiles
    given user input #x, list # mutated lines/gene'''

    '''
    Store the number of lines mutated for each gene,
    keyed by locus_tag
    '''
    tag_lines = {}

    for key, gd in genomediffs.items():
        if gd.mut_type == 'SNP' and gd.snp_type in snp_types:
            tag = tuple(gd.locus_tag)
            if tag not in tag_lines:
                tag_lines[tag] = {'genes': set(), 'lines': set(), 'gene_product': set()}
            tag_lines[tag]['lines'].add(gd.line)
            tag_lines[tag]['genes'].add(tuple(gd.gene_name))
            tag_lines[tag]['gene_product'].add(tuple(gd.gene_product))

    sorted_tag_lines = list(zip(list(tag_lines.keys()), list(tag_lines.values())))
    sorted_tag_lines = sorted(sorted_tag_lines,
                              key=lambda x: len(x[1]['lines']), reverse=True)

    return sorted_tag_lines


# TODO: Update for refactor? low priority
# TODO: if position is intergenic (or etc.), choose new position
def snpmutate(mat, num_replicates, refseq_arr):
    '''input: matrixdict and refseq as numpy array from
           snpcount and parse_ref, respectively
    called by get_mut_sites

    returns a dict with original base as key, numpy array as value:
    { 'A': np.array(num_replicates by num_mutatations A to others),
      'T': np.array(num_replicates by num_mutations T to others),
      ...
    }
    Where each array has num_replicates rows and columns
        corresponding to the number of mutations
    for the original base
    '''
# below was recently commented out for purposes of functionality
# mut_sites = {}
# for origbase in mat:
# num_muts = #numpy.zeros sum across matrix row

    '''
    mut_sites = {}
    #print matrix
    for origbase in matrix:
        num_muts = sum(matrix[origbase][newbase] for newbase
            in matrix[origbase])
        #print origbase, num_muts
        sites = np.zeros((num_replicates, num_muts), dtype=int)
        for rep in xrange(num_replicates):
            if rep % 50 == 0 and rep > 0: #prints progress report
                print '\t... rep', rep, 'for base', origbase
            try:
                rep_sites = np.random.choice(np.where(
                   refseq_arr==origbase)[0], size=num_muts, replace=False)
     # [0] because returns a tuple, and we just want 1st
        element (list of indices).
            except ValueError as e:
                print >>sys.stderr, e
                print >>sys.stderr, 'Error getting sites for replicate', rep
            else:
                sites[rep,:] = rep_sites
        mut_sites[origbase] = sites
    #print mut_sites
    return mut_sites
    '''


# TODO: Refactor? low priority. do for comparison btwn analytic and simuln.
def get_mut_sites(matrices, refseq, num_replicates):
    mut_sites = {}
    refseq_arr = np.array([c for c in refseq])
    for filename in matrices:
        print('\n** generating mutation sites for', filename)
        mut_sites[filename] = snpmutate(matrices[filename],
                                        num_replicates, refseq_arr)
    print()
    '''
    mut_sites = { 'filename1': { 'A': array(row for reps, columns
        for mut position indices),
                                 'T': array(...) ... },
                  'filename2': {...},
                  ...
                }
    '''
    return mut_sites


def SNPsToAlignment(conf):
    ''' rows are lexicographically sorted conf.GENOMEDIFF_FILES, and
    the last row is the reference sequence. columns are all positions
    that evolved in the set of genomes.
    '''

    ref_record = utils.parse_genbank(conf.REF_GENOME)
    # utils.print_genbank_summary(ref_record)

    snps = []
    # # each elt in snps is a tuple: (position, old_base,
    # new_base, locus_tag, label)
    conf.GENOMEDIFF_FILES.sort()  # # So I can assume the diffs are sorted.
    for gd_file in conf.GENOMEDIFF_FILES:
        gd_dict = parse_genomediff(gd_file, ref_record)
        for k, v in gd_dict.items():
            if v.mut_type != 'SNP':  # # only consider SNPs
                continue
   
         # old_base = ref_record[v.position]
            snps.append(
                (v.position + 1,
                 v.old_base,
                 v.new_base,
                 v.locus_tag[0],
                 gd_file))
    snps.sort(key=lambda elt: elt[0])  # sort by position.
    # cols = sorted([x for x in set([elt[0] for elt in snps])])
    ## NOTE: parse_genomediff converts 1-based indexing to 0-based indexing; 
    ## this line changes it back for reporting to be consistent with original gd files.
    cols = [x for x in set([elt[0] for elt in snps])]
    cols.sort()
    ## The LAST row of the alignment is the reference.
    alignment = [[''] * len(cols)
                  for x in range(len(conf.GENOMEDIFF_FILES)+1)]

    for elt in snps:
        i = conf.GENOMEDIFF_FILES.index(elt[4])
        j = cols.index(elt[0])
        alignment[-1][j] = elt[1]  # # the reference sequence.
        alignment[i][j] = elt[2]
    # now fill the empty entries in the matrix w/ the ref seq value.
    ref = alignment[-1]
    #print ref
    for i in range(len(alignment)):
        for j in range(len(cols)):
            if alignment[i][j] == '':
                alignment[i][j] = ref[j]

    str_alignment = [''.join(x) for x in alignment]
    aln_ids = [os.path.splitext(gd)[0] for gd in conf.GENOMEDIFF_FILES]
    aln_ids = aln_ids + [ref_record.id]  # add the reference.
    site_recs = [
        SeqRecord(
            Seq(x), id=y) for x, y in zip(
            str_alignment, aln_ids)]
    # # turn into a Biopython Alignment object.
    msa = MultipleSeqAlignment(site_recs)
    ## return both the msa as well as the position and gene for each column in the alignment.
    msa_annotation = []
    for i,pos in enumerate(cols):
        locus = None
        for elt in snps:
            if elt[0] == pos:
                locus = elt[3]
                break
        annotation = (i,pos,locus)
        msa_annotation.append(annotation)
    return msa, msa_annotation


def AlignmentToPhylogeny(aln):
    # # write alignment to a temporary file.
    out_handle = open("temp/aln.phy", "w")
    AlignIO.write(aln, out_handle, "phylip-relaxed")
    out_handle.close()
    # # raxml needs phylip formatted data.
    os.chdir("temp")  # # want RaxML output to go in the temp directory.
    raxml_cline = RaxmlCommandline(sequences="aln.phy",
                                   threads=2, model="GTRGAMMA", name="test")
    raxml_cline()
    os.chdir("..")
    tree = Phylo.read("temp/RAxML_bestTree.test", "newick")
    return tree


## from Biopython Phylo cookbook for looking up a Clade's parent.
def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child] = clade
    return parents

## hand-rolled postorder traversal order for the genotype tree (for Sankoff cost assignment)
## the code does not assume that we are dealing with a strictly binary tree.
def postorder_gtree(gtree,cur=0,stack=[]):
    ## if the cur node is not terminal, recur on children.
    if len(gtree[cur]['children']):
        stack_of_stacks = []
        for c in range(len(gtree[cur]['children'])):
            stack_of_stacks.append(postorder_gtree(gtree, gtree[cur]['children'][c],stack))
        return reduce(lambda x, y:x+y, stack_of_stacks) + [cur]
    else:
        return [cur]
    
## implementation of the Sankoff algorithm.
## assigns genotypes to internal nodes of the tree, and counts minimum mutations in tree.
def Sankoff(gtree):
    ## in future, might allow gap character state in alphabet, or other states
    ## representing different kinds of mutations. 
    alphabet = ['A','C','G','T']
    alnlen = len(gtree[postorder_gtree(gtree)[0]]['genotype'])
    for i in postorder_gtree(gtree):
        if not len(gtree[i]['children']): # is a leaf
            gtree[i]['cost'] = [{l:float("inf") for l in alphabet} for x in range(alnlen)]
            for index,letter in enumerate(gtree[i]['genotype']):
                gtree[i]['cost'][index][letter] = 0
        else: # is an internal node.
           gtree[i]['cost'] = [{l:0 for l in alphabet} for x in range(alnlen)]         
           for c in gtree[i]['children']:
               for x in range(alnlen):
                   child_cost_dict = gtree[c]['cost'][x]
                   best_child_state = min(child_cost_dict,key=child_cost_dict.get)
                   best_child_cost = child_cost_dict[best_child_state]
                   for l in alphabet: ## This is the fitch recurrence relation:
                       state_child_cost = child_cost_dict[l]
                       gtree[i]['cost'][x][l] += min(best_child_cost+1,state_child_cost)
    ## nodes are already in preorder, which is nice for backtracking.
    for i in gtree:
        if len(gtree[i]['genotype']) == 0: ## assign the genotype to the internal node.
            gtree[i]['genotype'] = ''.join([min(x,key=x.get) for x in gtree[i]['cost']])
    return gtree

def SankoffGenotypeTree(phy, aln, conf):
    ''' takes a Biopython Tree object and alignment. returns a tree of genotypes, with mutations
        stored in the internal nodes, using the Sankoff algorithm.
        This implementation uses a dict to represent a node, rather than a class.'''
    ## Walk down the Biopython Tree, and initialize the genotype tree.
    genotype_tree = {}
    clade_dict = {}
    parent = all_parents(phy)
    for i,clade in enumerate(phy.find_clades(order='preorder')):
        clade_dict[clade] = i
        if clade.is_terminal():
            ## assign genotype to the leaf.
            my_genotype = [rec for rec in aln if rec.id == clade.name][0].seq
            genotype_tree[i] = {'children':[], 'genotype':my_genotype, 'name':clade.name}
        else:
            genotype_tree[i] = {'children':[], 'genotype':''}
        if i > 0: # if not root, add current node as a child.
            parent_index = clade_dict[parent[clade]]
            genotype_tree[parent_index]['children'].append(i)
    ## Now, run the Sankoff algorithm on the tree with labeled leaves.
    return Sankoff(genotype_tree)


def BasicSNPCount(conf):
    ''' return the total number of nonsynonymous SNPs in genes,
        assumes all mutations are
        independent (star phylogeny).'''
    snpcount = 0
    ref_record = utils.parse_genbank(conf.REF_GENOME)
    for gd_file in conf.GENOMEDIFF_FILES:
        gd_dict = parse_genomediff(gd_file, ref_record)
        for k, v in gd_dict.items():
            if v.mut_type != 'SNP':  # # only consider SNPs,
                continue
            if v.snp_type == 'nonsynonymous':  # and those in genes.
                snpcount = snpcount + 1
    return snpcount


def formGenomeCDF(ref_record, loci):
    genes = []
    pdf = []
    all_genes_length = 0
    for f in ref_record.features:
        if f.type == 'CDS':
            all_genes_length = all_genes_length + len(f)
            locus_tag = f.qualifiers['locus_tag'][0]
            if locus_tag in loci:
                genes.append(locus_tag)
                pdf.append(len(f))
    pdf = [float(x) / float(all_genes_length) for x in pdf]
    cdf = []
    for i, p in enumerate(pdf):
        if i == 0:
            cdf.append(p)
        else:
            cdf.append(p + cdf[-1])
    return (genes, cdf)

        
def Statisticulate(conf, snptotal, tree_and_annotation=None, reps=1000):
    ''' Calculates statistics of parallel evolution given where mutations occurred. 
    right now, works only for nonsynonymous mutations. '''
    ref_record = utils.parse_genbank(conf.REF_GENOME)
    ######## count nonsynonymous mutations in each gene in the gd files.
    if tree_and_annotation is None: # default: assume star phylogeny.
        nonsynonymous_mutations = {}
        for gd_file in conf.GENOMEDIFF_FILES:
            gd_dict = parse_genomediff(gd_file, ref_record)
            for mut_id, gd in gd_dict.items():
                if gd.mut_type == 'SNP' and gd.snp_type == 'nonsynonymous':
                    # These all only have one locus_tag, and we can't
                    # use a list as a key, so just get the value
                    locus_tag = gd.locus_tag[0]
                    if locus_tag not in nonsynonymous_mutations:
                        nonsynonymous_mutations[locus_tag] = 1
                    else:
                        nonsynonymous_mutations[locus_tag] += 1

    else:  # a tree and annotation of mutation is provided.
        gtree, col_annotation = tree_and_annotation
        assert gtree is not None
        assert col_annotation is not None
        nonsynonymous_mutations = {}
        ## the root of gtree contains information about independent mutations.
        ## the cost at position X at the root is the number of independent mutations
        ## at position X (based on the given phylogeny).
        mut_counts = [min(x.values()) for x in gtree[0]['cost']]
        for i, mut_tuple in enumerate(col_annotation):
            column, pos, locus_tag = mut_tuple
            if locus_tag not in nonsynonymous_mutations:
                nonsynonymous_mutations[locus_tag] = mut_counts[i]
            else:
                nonsynonymous_mutations[locus_tag] += mut_counts[i]

    # # How many times does a dN occur in the gene?
    pval_numerator = {k: 0 for k in nonsynonymous_mutations}
    genes, cdf = formGenomeCDF(ref_record, list(nonsynonymous_mutations.keys()))
    for replicate in range(reps):
        for m in range(snptotal):
            nulldist = {k: 0 for k in nonsynonymous_mutations}
            # # draw a random number, and see which gene mutated.
            rando = random.random()
            if rando <= cdf[-1]:  # # rando is in the gene set.
                for i, x in enumerate(cdf):
                    if rando <= x:
                        nulldist[genes[i]] = nulldist[genes[i]] + 1
                        break  # # found the right bin.
        for g in nulldist:
            if nulldist[g] >= nonsynonymous_mutations[g]:
                pval_numerator[g] = pval_numerator[g] + 1
    pvals = {k: float(v) / float(reps) for k, v in pval_numerator.items()}
    for k, v in pvals.items():
        print("locus_tag:", k, "p-value:", v)

# def cleanup(): ## remove temp folder with inputs, and recreate it empty.
    # # THIS FUNCTION DOES NOT WORK RIGHT NOW
#    shutil.rmtree('./temp')
#    os.makedirs('temp')


def window_iterator(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(itertools.islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result


def makeWindow(i, mutlist, distlist, window_len=200):
    assert(len(mutlist) == len(distlist))
    x = i
    total_dist = 0
    window = []
    while total_dist <= window_len:
        window.append(mutlist[x])
        total_dist = total_dist + distlist[x]
        if x == len(mutlist) - 1:  # chromosome is circular.
            x = 0
        else:
            x = x + 1
    return window


def makeWindows(ref_record, mut_list, window_len=200):
    mut_list.sort(key=lambda x: x.position)
    dist_list = [y.position-x.position for x, y in window_iterator(mut_list)]

    # # add the last elt of dist_list (bc chromosome is circular)
    dist_list.append(mut_list[0].position +
        len(ref_record.seq) - mut_list[-1].position)
    assert(sum(dist_list) == len(ref_record.seq))
    windows = []
    for i, mut in enumerate(mut_list):
        windows.append(makeWindow(i, mut_list, dist_list, window_len))

    # # filter out structural mutations.
    windows2 = []
    for w in windows:
        kinds = set([x.mut_type for x in w])
        if kinds == set(['SNP']):
            windows2.append(w)
    windows2.sort(key=lambda x: len(x), reverse=True)

    return windows2


def pickWindows(conf, windows2):
    # # pick non-overlapping windows that cover all genomes as markers
    unmarked_genomes = {k:1 for k in conf.GENOMEDIFF_FILES}  # 1 if unmarked
    markers = []

    while sum(unmarked_genomes.values()):
        best_window = []
        most_new_marks = 0
        for w in windows2:
            # # make sure mutations UNIQUELY IDENTIFY genome.
            # # CANNOT be the same new_base and position.
            unique_positions = [mut.position for mut in w]
            id_check = [(mut.position, mut.new_base) for mut in w]
            unique_muts = [mut for mut in w if
                           id_check.count((mut.position, mut.new_base)) == 1]
            marks = [mut.fname for mut in unique_muts]
            new_marks = len([x for x in marks if unmarked_genomes[x] == 1])

            if new_marks > most_new_marks:
                most_new_marks = new_marks
                best_window = w

        markers.append(best_window)
        for mut in best_window:
            unmarked_genomes[mut.fname] = 0
                # the genome is not unmarked anymore!

    return markers


def printWindows(markers):
    # # print the ranges that are most informative.
    for i, m in enumerate(markers):
        print("MARKER", i+1, ":")
        for mut in m:
            print(mut.fname, mut.position, mut.old_base, \
                mut.new_base, mut.gene_name)


# TODO: Update for refactor
def write_gene_mut_counts(genecoords, mut_sites):
    header = 'gene, ' + ', '.join([filename for filename in mut_sites])
    with open('simulated_mutations_counts.csv', 'wb') as outfp:
        outfp.write(header + '\n')
        for locus_tag in genecoords:
            cds = genecoords[locus_tag]
            for filename in mut_sites:
                line_muts = 0
                for origbase in mut_sites[filename]:
                    line_muts += (
                        (mut_sites[filename][origbase] >= start) & (
                            mut_sites[filename][origbase] < end)).sum()
                row.append(line_muts)
            outfp.write(', '.join([str(c) for c in row]) + '\n')


# TODO: Update for refactor
def write_gd_gene_mut_counts(genecoords, gd_genes):
    header = 'gene, count'
    with open('experimental_mutations_counts.csv', 'wb') as outfp:
        outfp.write(header + '\n')
        for _, _, gene, tag in genecoords:
            muts = [tag]
            if tag in gd_genes:
                muts.append(gd_genes[tag])
            else:
                muts.append(0)
            outfp.write(', '.join([str(c) for c in muts]) + '\n')


# TODO: Update for refactor
def write_proc6_locus_mut_counts(linesmut):
    header = 'locus_tag; genomes'
    with open('locus_mut_counts.csv', 'wb') as outfp:
        outfp.write(header + '\n')
        for row in linesmut:
            locus = row[0]
            genomes = row[1]
            outfp.write('{}; '.format(locus))
            outfp.write(', '.join([str(g) for g in genomes]) + '\n')
