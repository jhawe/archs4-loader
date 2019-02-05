configfile: "./config.json"

# define subdir depending on cancer status

if(config["filter_cancer"] == "True"):
  DRESULTS = config["data_dir"] + "cancer_filtered/"
else:
  DRESULTS = config["data_dir"]

# ------------------------------------------------------------------------------
# Download the data file and extract the  global sample design sheet
# ------------------------------------------------------------------------------
rule init_data:
	output:
		"results/design_all_samples.tsv",
		"results/human_matrix_download.h5"
	conda:
		"envs/r_env.yaml"
	script:
		"scripts/create_sample_sheet.R"

# ------------------------------------------------------------------------------
# Extract and process expression data for specific keywords from main h5 file
# ------------------------------------------------------------------------------
rule extract_data:
	input:
		h5="results/human_matrix_download.h5"
	output:
		expr=DRESULTS + "{keywords}/expression_normalized.tsv",
		raw=DRESULTS + "{keywords}/expression_raw.tsv",
		design=DRESULTS + "{keywords}/design.tsv",
		plot=DRESULTS + "{keywords}/tsne.pdf"
	conda:
		"envs/r_env.yaml"
	wildcard_constraints:
		keywords="[a-zA-Z_]+"
	params:
		norm_method = "sva", # must be one of 'sva' or 'quantile'
		filter_cancer = config["filter_cancer"]
	log:
		"logs/extract_data/{keywords}.log"
	script:
		"scripts/extract_data.R"

# ------------------------------------------------------------------------------
# Extract data from ARCHS4 with exact keyword matching
# We use this to be able to define certain subsets of samples for further 
# processing, e.g. based on intitial fuzzy matching results
# ------------------------------------------------------------------------------
tissue_to_keyword = {"pancreas":"pancreas|human pancreas",
  "liver":"Human liver tissue|liver biopsy|Liver|normal liver|Normal liver|normal liver tissue|liver tissue|liver",
  "whole_blood":"Whole Blood|whole blood|Whole blood|whole peripheral blood|whole blood cells|peripheral whole blood",
  "lung":"lung"}


rule extract_data_exact:
	input:
		h5="results/human_matrix_download.h5"
	output:
		expr=DRESULTS + "exact/{tissue}/expression_normalized.tsv",
		raw=DRESULTS + "exact/{tissue}/expression_raw.tsv",
		design=DRESULTS + "exact/{tissue}/design.tsv",
		plot=DRESULTS + "exact/{tissue}/tsne.pdf"
	conda:
		"envs/r_env.yaml"
	wildcard_constraints:
		keywords="[a-zA-Z_]+"
	params:
		keywords=lambda wildcards: tissue_to_keyword[wildcards.tissue],
		filter_cancer = config["filter_cancer"]
	log:
		"logs/extract_data_exact/{tissue}.log"
	script:
		"scripts/extract_data_exact.R"

# ------------------------------------------------------------------------------
# Explore the data and create a nice summary
# Note: only single output files can be used for R markdowns
# ------------------------------------------------------------------------------
rule explore_data:
	input:
		expr=DRESULTS + "{keywords}/expression_normalized.tsv",
		raw=DRESULTS + "{keywords}/expression_raw.tsv",
		design=DRESULTS + "{keywords}/design.tsv"
	output:
		DRESULTS + "{keywords}/summary.html"
	params:
		pdf_out = DRESULTS + "{keywords}/tsne_plots.pdf"
	conda:
		"envs/r_env.yaml"
	script:
		"scripts/explore_data.Rmd"

# ------------------------------------------------------------------------------
# Target rule to process some interesting tissues
# ------------------------------------------------------------------------------
def gtex_tissue_files(wc):
	with open("gtex_tissues.txt") as f:
		tissues = f.readlines()
	tissues = [x.strip() for x in tissues]
	return(expand(DRESULTS + "{keywords}/summary.html", keywords=tissues))

rule get_gtex_tissue_summaries:
	input:
		gtex_tissue_files
	output:
		"results/summaries.zip"
	shell:
		"zip {output[0]} {input}"

rule explore_gtex:
	input:
		tissues="gtex_tissues.txt",
		h5="results/human_matrix_download.h5"
	output:
		expr=DRESULTS + "all_gtex/expresion_normalized.tsv",
		raw=DRESULTS + "all_gtex/expresion_raw.tsv",
		design=DRESULTS + "all_gtex/design.tsv",
		plot=DRESULTS + "all_gtex/tsne.pdf"
	threads: 10
	params:
		norm_method = "sva"
	log:
		"logs/explore_gtex.log"
	script:
		"scripts/explore_gtex.R"
