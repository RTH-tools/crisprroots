##Version 1.1
### Technical improvements
  * Substitution of CRISPRoff v.1.1.1 to CRISPRoff v.1.1.2, ported to Python3.
  * Removal of environment py2, not necessary for CRISPRoff v.1.2.
  * Installation of the ViennaRNA package moved from py2 to py3.
  * Removal of R packages installations from R environment. Now the requirements for the standard
  workflow are satisfied in Conda. 
  * Automatization of the path setting in the confing file of the test dataset
### Bugs fixed
  * Environment twobit was failing if Python3 was not preinstalled. Installation of Python3 is now added to environment twobit.
  * Fixed a bug in the rule "GATK_mutect2_chromosome_wise" which replaced "." in folder names with "_"

### Other changes
  * Cleaning of the config file (comments, unused variables)
  * Changes in the variables names in the config
  * Changes to the test datasets: weight for GG PAM binding site is 1.0, not 0.1
  * Changes to the test resources: lighter STAR indexed chromosome (only the portion of the chromosome to which reads map is indexed)
  * The sheet name in the on-target_knockin.xlsx output is changed from "Detailed" to "Summary KI"
  * Additional output files are marked as temp. See the Manual for a full overview of which output files are temp and which are always output by the pipeline
  * Addition of the folder replicate_publication_results with the scripts necessary to replicate all of the results presented in the CRISPRroots publication
