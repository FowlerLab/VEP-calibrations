# Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes

  This github repo contains code to re-create the analysis and figures from "Calibration of variant effect predictors on genome-wide data masks heterogeneous performance across genes". 
  Supplementary files from the paper can be found at https://zenodo.org/records/11256843

## Analysis and Figures 

  The analysis and figures folder contains scripts used for the paper's analysis and generating the figures in the paper. In most cases, you will need the supplemental files to run the script, in some instances, files are provided with the script. The beginning of each script describes the files you will need to run it. 

  ### Figure 1 plots
  This folder contains the script to generate the whole genome histograms seen in Figure 1. It is not recommended to run this script locally as the supplemental files needed for the script are very big.

 ### Figure 2, Figure 5, and Supplemental Figures 1 and 2 and gene-specific analysis 

   This folder contains two scripts, one for the REVEL predictor and one for the BayesDel predictor. Each script performs the gene-specific analysis of concordance/discordance from the paper and generates a PDF containing Figure 2 (or supplemental Figure 1) and Figure 5 (or supplemental Figure 3). These scripts create multiple intermediate data frames with specific analyses; any data frames can be downloaded for further analysis using the write.csv() command. 

### Figure 3

  This folder contains the script to create Figure 3 from the paper and any files needed to run the script are on Zenodo (linked above)

### Figure 4 

  This folder contains the script to create Figure 4 histograms from the paper, there are two separate sections of code to produce a REVEL histogram and a BayesDel histogram. The necessary files to run the script are in the folder. The script creates a pdf of the six plots from Figure 4. 

## Parse

The parse folder contains scripts that filter the original ClinVar download annotated with predictor scores using Annovar. Scripts produce intermediate .csv files for further filtering. The final script filtering scripts create the ClinVar 2023 Dataset, the ClinVar 2023 dataset without training variants, and the ClinVar VUS dataset. The original annovar annotate ClinVar file with predictor scores can be found here: https://zenodo.org/records/11256843  

  ### Annovar parse
  
  This folder contains a Python script to parse through the ClinVar annotated file with predictor scores from Annovar. It will create intermediate files that are nice to have but not necessary. Only two files needed for downstream filtering will be produced from dfa (for Clinvar 2023 VUS dataset) and df6 (Clinvar 2023 dataset and Clinvar 2023 dataset without training variants)

### Filtering 

  This folder contains an R script to filter for allele frequency, remove gene duplicates, and filter for disease-relevant genes to create the Clinvar 2023 dataset. It also uses the file from dfa (earlier) and filters it to remove gene duplicates to create the ClinVar 2023 VUS dataset. Additional relevant files will also be found in the folder. 

### Remove training variants 

This folder contains an R script to filter out the training variants from the Clinvar 2023 dataset to create the Clinvar 2023 dataset without training variants. Additional relevant files will be found in the folder. 


   
