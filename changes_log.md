##Version 1.3
### Technical improvements

### Bugs fixed
  * The input-output connection between steps 5- and 7- was incorrect, and because of that all samples had to run together, in a single job, in step 5. Now they are correctly parallelized.
  * The expression-based off-target screening is modified to remove potential off-targets that cleave the on-target site. In "find_DE_genes" the on-target KO gene is not removed from the list of DE genes that can carry off-targets.
  * Corrected wrong table names in get_offtarget_context_for_sequencing.py
### Other changes
  * In the variant-based off-target screening, the "lifted edits" file is opened before the loop over variants, to avoid multiple opening operations
  * Threshold for DE genes changed to 0.05 instead of 0.01 in expression-based off-target search
  * Introduced new columns in report files to state the full gRNA-pOT binding pattern and allele counts that support variants linked tovariant-based pOTs (pOT=potential off-target)
  * Terminated eSNPKaryotyping support
### Environment changes
  * Terminated conda support and switched to docker/singularity for improved stability. The definition of the conda environments is still present and functional. However, due to the instability of certain environments, we will no longer offer support on this and recomment using the provided docker container.

