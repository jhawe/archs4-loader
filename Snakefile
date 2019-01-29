configfile: "./config.json"

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
		expr=config["data_dir"] + "{keywords}/expression_normalized.tsv",
		raw=config["data_dir"] + "{keywords}/expression_raw.tsv",
		design=config["data_dir"] + "{keywords}/design.tsv"
	log:
		"logs/extract_data_{keywords}.log"
	conda:
		"envs/r_env.yaml"
	params:
		norm_method = "sva" # must be one of 'sva', 'peer' or 'quantile'
	log:
		"logs/extract_data/{keywords}.log"
	script:
		"scripts/extract_data.R"

# ------------------------------------------------------------------------------
# Explore the data and create a nice summary
# Note: only single output files can be used for R markdowns
# ------------------------------------------------------------------------------
rule explore_data:
	input:
		expr=config["data_dir"] + "{keywords}/expression_normalized.tsv",
		raw=config["data_dir"] + "{keywords}/expression_raw.tsv",
		design=config["data_dir"] + "{keywords}/design.tsv"
	output:
		config["data_dir"] + "{keywords}/summary.html"
	params:
		pdf_out = config["data_dir"] + "{keywords}/tsne_plots.pdf"
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
	return(expand(config["data_dir"] + "{keywords}/summary.html", keywords=tissues))

rule explore_gtex:
	input:
		gtex_tissue_files
	output:
		"results/summaries.zip"
	shell:
		"zip {output[0]} {input}"
