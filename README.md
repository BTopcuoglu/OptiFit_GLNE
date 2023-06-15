## Machine learning classification by fitting amplicon sequences to existing OTUs

The ability to use 16S rRNA gene sequence data to train machine learning classification models offers the opportunity diagnose patients based on the composition of their microbiome. In some applications the taxonomic resolution that provides the best models may require the use of de novo OTUs whose composition changes when new data are added. We previously developed a new reference-based approach, OptiFit, that fits new sequence data to existing de novo OTUs without changing the composition of the original OTUs. While OptiFit produces OTUs that are as high quality as de novo OTUs, it is unclear whether this method for fitting new sequence data into existing OTUs will impact the performance of classification models relative to models trained and tested only using de novo OTUs. We used OptiFit to cluster sequences into existing OTUs and evaluated model performance in classifying a dataset containing samples from patients with and without colonic screen relevant neoplasia (SRN). We compared the performance of this model to standard methods including de novo and database-reference-based clustering. We found that using OptiFit performed as well or better in classifying SRNs. OptiFit can streamline the process of classifying new samples by avoiding the need to retrain models using reclustered sequences.

### Citation
> Armour CR, Sovacool KL, Close WL, Topçuoğlu BD, Wiens J, Schloss PD.
> 2023. Machine learning classification by fitting amplicon sequences to existing OTUs
> (URL)
### Quickstart
1. Clone this repository.
   ```
   git clone https://github.com/SchlossLab/Armour_OptiFitGLNE_XXXX_2023
   ```
2. Install the dependencies.
   ```
   conda env create -f envs/glne.yaml
   conda activate glne
   ```
3. run the entire pipeline.  
   locally:
   ```
   snakemake --cores 4
   ```
   on an HPC running slurm:  
   (You will first need to edit your email and slurm account info in the
    [submission script](code/slurm/snakemake.sh)
    and [cluster config](config/config.yaml)
   ```
   sbatch code/slurm/snakemake.sh
   ```
   

### Directory Structure

	|- README          # the top level description of content (this doc)
	|- CONTRIBUTING    # instructions for how to contribute to your project
	|- LICENSE         # the license for this project
	|
	|- code/          		# any programmatic code
	| |- bash/		  		# bash scripts, including all mothur commands
	| |- R/			  		# R scripts
	| |- py/		  		# python scrips
	| +- slurm/snakemake.sh # script for submitting jobs to the cluster
	|
	|- config/
	| +-config.yaml   # cluster configuration file
	|
	|- data           # raw and primary data, are not changed once created
	| |- metadata/    # study metadata
	| |- process/	  # cleaned data
	| |- raw/         # raw data
	| +- references/  # reference files to be used in analysis
	|
	|- docs/
	| +- exploratory.html # document with additional exploratory analysis
	|
	|- envs/
	| +- glne.yaml    # conda environment file
	|
	|- exploratory/   # exploratory data analysis for study
	| |- notebook/    # preliminary analyses
	| +- scratch/     # temporary files that can be safely deleted or lost
	|
	|- results        # all output from workflows and analyses
	| |- figures/     # graphs, likely designated for manuscript figures
	| |- ml/	  	  # output from ML models
		|- results/	  # results from models for each split
		|- summary/   # merged results from all splits of models
	| +- tables/      # text version of tables to be rendered with kable in R
	|
	|- submission/
	| |- study.Rmd      # executable Rmarkdown for this study, if applicable
	| |- study.md       # Markdown (GitHub) version of the *.Rmd file
	| |- study.tex      # TeX version of *.Rmd file
	| |- study.pdf      # PDF version of *.Rmd file
	| |- header.tex     # LaTeX header file to format pdf version of manuscript
	| |- references.bib # BibTeX formatted references
	| |- XXXX.csl       # csl file to format references for journal XXX
	| +- figures/		# manuscript figures
	|
	+- Snakefile        # executable Snakefile for this study
