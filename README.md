## Streamlined implementation of a machine learning model to classify screen relevant neoplasia using reference-based OTU clustering

Machine learning classification of disease based on the gut microbiome often relies on clustering 16S rRNA gene sequences into operational taxonomic units (OTUs) to quantify microbial composition. The abundance of each OTU is then used to train a classification model. The standard de novo approach to clustering sequences into OTUs leverages the similarity of the sequences to each other rather than to a reference database. However, such an approach depends on the sequences in the dataset and therefore OTU assignments can change if new data are added. This lack of stability complicates classification because in order to use the model to classify additional samples, the new sequences must be reclustered with the old data and the model must be retrained with the new OTU assignments. The new reference-based clustering algorithm, OptiFit, addresses this issue by fitting new sequences into existing OTUs. While OptiFit can produce high quality OTU clusters, it is unclear whether this method for fitting new sequence data into existing OTUs will impact the performance of classification models. We used OptiFit to cluster additional data into existing OTU clusters and then evaluated model performance in classifying a dataset containing samples from patients with and without colonic screen relevant neoplasia (SRN). We compared the performance of this model to the standard procedure of de novo clustering all the data together. We found that both approaches performed equally well in classifying SRNs. Moving forward, when OTUs are used in classification problems, OptiFit can streamline the process of classifying new samples by avoiding the need to retrain models using reclustered sequences.

### Citation
> Armour CR, Sovacool KL, Close WL, Topçuoğlu BD, Wiens J, Schloss PD.
> 2022. Streamlined implementation of a machine learning model to classify screen relevant neoplasia using reference-based OTU clustering. 
> https://doi.org/10.1101/2022.09.01.506299
### Quickstart
1. Clone this repository.
   ```
   git clone https://github.com/SchlossLab/Armour_OptiFitGLNE_mBio_2023
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

	project
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
