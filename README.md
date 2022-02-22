# CRISPR-Cas9-fitness-effects
Files and data to accompany "Fitness effects of CRISPR endonucleases in *Drosophila melanogaster* populations"

**AML-ML-04092020-emp-v10.R**: R-code to quantify the fitness costs of the four different transgenetic constructs from the observed frequency trajectories in our *D. melanogaster* cage populations. Uses functions implemented in AML-ML-v10.R.

**AML-ML-09042021-sim-v10.R**: R-code to simulate construct frequencies given different fitness cost estimates. Uses functions implemented in AML-ML-v10.R.

**AML-ML-v10.R**: R-code of previously developed maximum likelihood framework (*https://doi.org/10.1534/genetics.118.301893*) modified to model two unlinked autosomal loci, representing the construct and a single idealized off-target site. 

**gene_drive_off_target_effects.slim**: Agent-based simulation framework implemented in SLiM3 to study the consequences of off-target fitness costs on the invasion potential of homing and suppression gene drives. 

**constructs.zip**: annotated sequences of the final construct insertions in ApE format (*http://biologylabs.utah.edu/jorgensen/wayned/ape*).

**data/**: Data used by AML-ML-v10.R. FitnessModes.txt contains the tested models encoded as integer vectors; rawData/ contains the raw genotype counts of our *D. melanogaster* cage populations (can also be found in Supplemental_Data_Sets.xlsx); TransitionMatrix.txt contains expected offspring counts for each genotype combination.

**example-pictures.zip**: A picture sample set for the image-based screening pipeline (Macro-Feb2020-manualClassification.ijm).

**Macro-Feb2020-manualClassification.ijm**: Macro of the image-based screening pipeline. 

**phenotype_assays.zip**: Raw data of the phenotypic assays. Used by Phenotypes.R.

**Phenotypes.R**: R-code to analyze the phenotypic assays.

**Supplemental_Data_Sets.xlsx**: Raw counts of each experimental population (different constructs and the Cas9HF1 homing drive). 
