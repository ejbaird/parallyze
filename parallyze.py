#usage: python parallyze.py

import argparse
import os
import sys

import config
import config_default

def file_list(fs):
    flist = fs.split()
    for f in flist:
        assert os.path.isfile(f)
    return flist

def get_config():
    conf={}
    if config.PROCEDURE == '3':
        ref = config.REF_GENOME.strip()
        assert os.path.isfile(ref)
        assert ref.endswith('.gb') 
        conf['procedure'] = config.PROCEDURE
	conf['ref'] = ref
	assert len(config.GENOME_DIFFS.strip())==0
        return conf
    elif config.PROCEDURE in ['1','2','4','5']:  
        ref = config.REF_GENOME.strip()
        assert os.path.isfile(ref)
        assert ref.endswith('.gb')
        diffs = file_list(config.GENOME_DIFFS)
	assert len(diffs)!=0
        for diff_file in diffs:   #annotated vs. non-annotated genomediff files
	     assert diff_file.endswith('.gd')
        conf['ref'] = ref
        conf['procedure'] = config.PROCEDURE
        conf['diffs'] = diffs
        return conf
    else:
        print >>sys.stderr, 'Invalid procedure {p}! Using default data'.format(p=config.PROCEDURE)
        conf['ref'] = config_default.REF_GENOME
        conf['procedure'] = config_default.PROCEDURE
        conf['diffs'] = config_default.GENOME_DIFFS
        return conf
    print class (conf['diffs'])

get_config()

#fp=open(filename, 'rU')

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq

from numpy import random

def base_to_int(i):
    if i == 'A':
        return 0
    elif i == 'G':
        return 1
    elif i == 'C':
        return 2
    else:
        return 3

def int_to_base(i):
    if i == 0:
        return 'A'
    elif i == 1:
        return 'G'
    elif i == 2:
        return 'C'
    else:
        return 'T'

def seq_to_int(seq):
    converted=[base_to_int(b) for b in seq]
    return converted

def int_to_seq(seq):
    converted=[int_to_base(b) for b in seq]
    return converted

def proc1(conf):
    mutlist=[]
    mutations={}
    with open(conf['diffs']) as fp: #conf['diffs'] should be string or buffer?
        for line in fp: 
            if line.startswith('#' or 'JC' or 'RA' or 'UN'):
                continue #does this start 'for line' loop again?
            line=line.split()
            mut_type=line[0]
            mut_id=line[1]
            parent_ids=line[2].split(',')
            seq_id==line[3]
            position==line[4]
            if mut_type=='SNP':
                new_seq==line[5]
            elif mut_type=='SUB':
                size==line[5]
                new_seq==line[6]
            elif mut_type=='DEL':
                size==line[5]
            elif mut_type=='INS':
                new_seq==line[5]
            elif mut_type=='MOB':
                repeat_name==line[5]
                strand==line[6]
                duplication_size==line[7]
            elif mut_type=='AMP':
                size==line[5]
                new_copy_number==line[6]
            elif mut_type=='CON':
                size==line[5]
                region==line[6]
            elif mut_type=='INV':
                size==line[5]
            #while mut_type count is 3 chars, continue - to terminate transfer of info when it gets gobbledy-gooky?

    #mutations[mut_id]={'type': mut_type, 'parents': parent_ids, ...)
    #what is funciton of above line?
def proc3(conf):
    print '\n', 'Assumptions:', '\n', 'Synonymous mutations are neutral' '\n', 'Infinite sites model', '\n', 'Mutations are independent of one another', '\n', 'No defects to DNA repair', '\n', 'Mutation rate is constant across the genome', '\n', 'There is only one chromosome', '\n'
    lines=input("How many lines?  ")
    gens=input("How many generations? ")
    reps=input("How many replicates?  ")

    '''
    SeqIO.parse is an iterator and so has a method named "next"
    which is called when you use it in a for loop, ie:
    for record in SeqIO.parse(stuff):
        do stuff
    We can call it explicitly, once, because we're assuming there is only
    one record in the SeqIO iterator. We will be *very* explicit and store
    the iterator itself as it, then  call next() on it like so:
    '''
    it = SeqIO.parse(conf['ref'], "genbank")
    record = it.next() 

    # convert the biopython Seq object to a python string
    seq = list(str(record.seq))
    length=len(seq)
    ## print 'Seq as list [truncated]:', seq[:5000], '...'
    
    countA=seq.count('A') #possibly change to 0,1,2,3
    countG=seq.count('G')
    countC=seq.count('C')
    countT=seq.count('T')
    print 'A:', countA, '   G:', countG, '   C:', countC, '   T:', countT

    # or do seq = ''.join(seq) to save it as a string and overwrite the list
    ##print 'Seq as string [truncated]:', ''.join(seq)[:1000], '...'
    print "Number of bases: ", length
    print 'Seq as condensed string:', ''.join(seq)[0:100], '...', ''.join(seq)[length-100:length] 

def main():
    parser = argparse.ArgumentParser()
    #parser.add_argument(dest='config', help="Config.py file should be in working folder.")
    #parser.add_argument('-r', dest='reference', default="NONE")
    args = parser.parse_args()

    conf = get_config()
    if conf['procedure']=='3':
	print >>sys.stderr, 'Configuration', '\n', 'Procedure: ', conf['procedure'], '\n','Reference: ', conf['ref']
	proc3(conf)
    else: 
        print >>sys.stderr, 'Configuration', '\n', 'Procedure: ', conf['procedure'], '\n','Reference: ', conf['ref'], '\n','Genome diffs: ', conf['diffs']
	if conf['procedure']=='1':
	    proc1(conf)
	elif conf['procedure']=='2':
	    proc2(conf)
	elif conf['procedure']=='4':
	    proc4(conf)
	elif conf['procedure']=='5':
	    proc5(conf)
main()
