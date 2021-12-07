##Version 1.2
### Technical improvements

### Bugs fixed
  * The sorting of potential off-targets based first on RISK, then on DeltaG_B condisered DeltaG_B as a string, rather than a float number. 
  * Variants overlapping deletions (*) were integrated in the genome during the generation of the variant-aware reference, but were automatically removed by faToTwoBit, resulting in a shift of coordinates between the fasta and the two-bit variant-aware genomes. This is now solved by transforming \* to N. This does not affect the results presented in the CRISPRroots publication.  

### Other changes
  * New dependency for jinja2 in the py3 conda environment, required for styling the candidate off-targets output file
  * The output table of candidate off-targets has been improved with colors highlighting specific entries. The excel file now include a front page with information on how to interpret the content of the tables.
  * New script on how to obtain DNA sequences flanking potential off-targets, which can be used to design primers. The script is now described in the manual.

