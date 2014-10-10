parallyze
=========

You will need to install Biopython before running this program. Either do the following from the command line or download BioPython directly from the website: http://biopython.org/wiki/Download 
        
    sudo pip install biopython
    sudo pip install numpy

Before running this program, ensure you have appropriately filled out your configuration file with a Genbank-formatted reference genome and annotated genome diff files (http://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing, 'ANNOTATE' command in gdtools). 

To run the program, type this into the command line:

    python parallyze.py

software to analyze parallel genome evolution by generating null distributions.

parallyze can also be used to do a power analysis for the number of lines 
needed in an evolution experiment to look for parallelism, given the mutation rate and spectrum, 
as well as the length of the experiment.

## Algorithm design

In all cases, the user should be able to specify only examining coding regions, or just nonsynonymous and synonymous mutations, or just synonymous mutations.

###Procedure 1: Gene mutation frequency. For each gene, what is the frequency of mutation for all of my pooled experimental lineages compared to a simulated null distribution? 

Genomes from evolution experiment. Assume all independent, i.e. star phylogeny.

Count all point mutations (x_1 + x_2 + ... + x_n), and turn into 4x4 matrix.

"Sort" reference genome positions by base. represent as (base, position), e.g. ('A', 3461000)

    for 1 to N replicates:
        for 1 to n genomes:
          draw x_i mutations from the mutation matrix, and drop onto reference genome.

###Procedure 2: Dispersion Test. For a single gene, what is the distribution of mutations across my experimental lineages compared to a null distribution? 

Genomes from evolution experiment. Assume star phylogeny.

    for 1 to N replicates:
        Shuffle all mutations across n genomes.
    Calculate how often a certain dispersion pattern occurs 
    (e.g., 12 mutations in nadR; all mutations occur in separate genomes).    
        
###Procedure 3: Phylogeny construction. What is the phylogeny of my sequenced genomes?

Genomes from multiple isolates from the same experimental evolution population.

1) Infer phylogeny
2) Infer genotypes of all internal nodes by "using parsimonious assumptions" -- or better.
3) Count the number x of inferred mutations, and generate a 4x4 matrix of mutation probabilities.

        for 1 to N replicates:
            drop x mutations onto reference genome, and count number of independent mutations per gene.
        
        average the results to calculate the null distribution.

This procedure could be extended to clinical or epidemiological isolates in the future.

###Procedure 4: What are the relative counts of dN, dS, and intergenic mutations at gene and genome level?

This should be straightforward from the genome diff format.

###Procedure 5: Gene mutation frequency. For each gene, what is the frequency of mutation for all of my pooled experimental lineages compared to an analytical null distribution? 

###Procedure 6: Number of lines mutating in each particular gene. /done - integrated into other procedures

###Datasets for Testing

* 40K Clones from LTEE
* The matrix of all LTEE lines sequenced over time.
* Bennett temperature-evolved genomes (doi:10.1126/science.1212986)
* Brian Wade's dessication lines
* Josh Nahum and Christian's sequencing of Paco's lines.
* Phage lambda datasets?
* Perhaps Lieberman et al. Burkholderia outbreak isolates (doi:10.1038/ng.997)?

##NOTES

The 12 lines of the LTEE probably don't have enough statistical power to
search for compensatory adaptation; but this software might be useful for designing experiments to detect
statistical signatures of compensatory adaptation. This depends on the native expectation of multiple hits in genes,
without consideration of protein stability.

In hypermutator lineages, genes with either zero OR multiple mutations might be candidates for compensation:
an underdispersion signal, as opposed to the overdispersion signal of strong parallelism.

Different mutational processes can be superimposed on each other (indels, rearrangements, transpositions).
Future extensions could code more sophisticated mutational models, trained on actual data from evolution experiments.

Look at Tenaillon paper for hierarchies of parallelism (gene level, operon level, pathway level, etc.)

Rich also mentioned tests on comparing the temporal order of mutations--to see if pykF always happens before spoT, for example. However, since it's difficult to calculate or know about the true target size (say if spoT is 5 times as large
compared to pykF, but if only 2 sites in the spoT matter, compared to pykF, etc.), this analysis may not fly.

In any case, including a temporal dimension (randomizing the identity of mutations over a phylogeny)
will surely allow for other interesting statistical tests.
