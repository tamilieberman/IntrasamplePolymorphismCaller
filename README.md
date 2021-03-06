Intrasample Polymorphism Caller
================


Some MATLAB functions to aid in finding polymorphic nucleotide positions from single-chromosome small genomes, not optimized for speed. <br> <br>
*pileup_to_diversity_matrix.m* summarizes pileup formats in a convenient 39 x genomelength matrix for easy interrogation in MATLAB (pileup_to_diversity_matrix.m) <br>
*find_diverse_positions_single_sample.m* deploys the filters in Lieberman et al, 2014 to call polymorphic positions.

It is strongly recommended that you use the interactive enviornment facilitated by MATLAB to compare these 39 attributes across your samples at each position. <br>

If you find these scripts helpful, please cite: <br>
Lieberman, T. D., Flett, K. B., Yelin, I., Martin, T. R., McAdam, A. J., Priebe, G. P., & Kishony, R. (2014). Genetic variation of a bacterial pathogen within individuals with cystic fibrosis provides a record of selective pressures. Nature Genetics.



Steps to analyzing intrasample diversity: <br>
------------------------------------------------------------
1) Summarize SAMtools pileup files in a convenient matlab structure using *pileup_to_diversity_matrix.m*  <br>
2) Find major allele frequencies in negative control sample using *maf_from_matrix.m*.  (optional) <br>
3) Find polymorphic positions in your sample using *find_diverse_positions_single_sample.m* <br>
4) Adjust filters by comparing across samples.



(1) Summarize SAMtools pileup files in a convenient matlab structure using *pileup_to_diversity_matrix.m*: <br>
------------------------------------------------------------
See header of this file for use on this command and description of the resulting data matrix.<br>
For a single 4 Mbp sample with 100x coverage, this function can take over an hour to run on a fast laptop --  we recommend running this on each of your files seperately on a cluster. It is slow largely because it peforms a huge number of statiscal tests, but can surely be optimized to run faster. <br>
Usage:<br>
pileup_to_diversity_matrix(INFILE, OUTFILE, GENOMELENGTH)<br>
Example:<br>
pileup_to_diversity_matrix('sample.pileup', 'sample.mat', 4411532) <br>

(2) Find major allele frequencies in negative control sample using *maf_from_matrix.m*: (optional)<br>
------------------------------------------------------------
One of the most useful filters is masking polymorphic positions found in a negative (isogenic) control sample. This function quickly calculates the major allele frequency at each genomic position in your negative control (or any sample) using the output from pileup_to_diversity.mat<br>
Usage:<br>
*control_freqs=maf_from_matrix(INFILE,OUTFILE)*<br>
Example:<br>
*control_freqs=maf_from_matrix('control.mat')*<br>


(3) Find polymorphic positions in your sample using *find_diverse_positions_single_sample.m*:  <br>
------------------------------------------------------------
Returns a numeric list of positions on the genome with polymorphisms. <br>
And example PARAMETERSFILE is included in this repository. Thresholds are all strictly greater than or less than, and should be self explanatory when combined with header of *pileup_to_diversity_matrix.m*.<br>

Usage:<br>
*positions=find_diverse_positions_single_sample(INFILE,PARAMETERSFILE,OUTPUTSUMMARYFILE,CONTROLFREQS)*<br>
or (if no isogenic control) <br>
*positions=find_diverse_positions_single_sample(INFILE,PARAMETERSFILE,OUTPUTSUMMARYFILE)*<br>
Example:<br>
*positions=find_diverse_positions_single_sample('sample.mat','example_params.mat','summary.txt',control_freqs)*<br>


(4) Adjust filters by comparing across samples. <br>
------------------------------------------------------------
Please write your own scripts in MATLAB, investigate the data carefully, and make sure that filters are set appropriately for your genome, samples, and application.<br>
It is recommended to grab all 39 attributes from the .mat files at each candidate variant position in your samples, to make a 3 dimensional matrix of dimension (39 X numpositions X numsamples). Filters can be set loosely at first, and then made more strict upon invesitgation of the data across samples.
