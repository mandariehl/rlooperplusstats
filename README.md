# rlooperplusstats
Code from Vasquez, Jonoska, et.al. paper "Modeling RNA:DNA Hybrids with Formal Grammars"

The file main.py:

Given an input string of nucleotides (DNA template strand read in 5' to 3' direction), we would like to check if the nearby region of nucleotides after are favorable for R loop formation. 


Output: binary string, 0 for not double C-rich (meaning not favorable for an r-loop) or 1 for double C-rich. We then use another function to determine if we have "enough" double-G-rich values in the neighborhood to call the region favorable for R-loop (1) or unfavorable (0). So it is a two-step process: First decide if each substring is double_C_rich, then decide if we have enough double_c_rich regions close enough together to start an R-loop.

## There are many helper functions that appear as part of the first-step function (which is double_c_rich), and several auxiliary functions that are used for determining the appropriate thresholds to use as parameters in double_c_rich.  

##double_c_rich(input_string, m, j, threshold1, threshold2)
# input_string: should be a template DNA sequence read 5' to 3'
# m: how often we would like to check for double_c_richness. m = 1 means we look at a substring starting at every nucleotide. m=5 checks only every fifth nucleotide. 
# j: how long a substring we would like to examine for double_c_richness. j = 10 means a moving window of size 10. 
#threshold1: this is the threshold in order to consider a substring as C-skewed. The calculation for C-skewness is (# of C's - # of G's)/(# of C's + # of G's), so a threshold1 of 0.3 means the C's to G's is at least 65% to 35%, or essentially a 2 to 1 ratio. 
#threshold2: this is the threshold in order to consider a substring as C-rich. This calculation is # of C's / # of nt. Note that C-skewness is not a purely stronger measure, since it doesn't include information about the proportion of A's and T's. A threshold of 0.6 means that at least 60% of the j nucleotides are C's. 

##The second-step function is findjonesoutoften(input_string, j).
# input-string: a binary string indicating the double_c_richness of a template DNA strand
# j: how many times we would need to see a 1 in a substring of length 10 before we consider it R-loop favorable (1) or unfavorable (0). 

#I could certainly combine those into a single function, but I prefer the two-step process at this point since we might still tune our thresholds as we get more data. 

The last functions in main.py are used for determining which thresholds match the output of the r-loop prediction software R-looper (written by Stolz as a member of the Chedin lab). We use spearmanr coefficient to tune our thresholds on two sample plasmid strands. The template strands in 5' to 3' direction are provided WITHIN the code as pfc53 and pfc 19. These strands were obtained via private communication fromt the Chedin lab. 


