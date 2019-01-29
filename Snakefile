configfile: "./config.json"

# ------------------------------------------------------------------------------
# Extract the  global sample design sheet from the h5 archive
# ------------------------------------------------------------------------------
rule create_sample_sheet:
	output:
		"results/design_all_samples.tsv",
		"results/human_matrix_download.h5"
	script:
		"scripts/create_sample_sheet.R"

# ------------------------------------------------------------------------------
# Get samples which have specific keywords annotated in the tissue-metadata
# ------------------------------------------------------------------------------
rule get_samples:
	input:
		design="results/design_all_samples.tsv",
	output:
		"results/sample_definitions/{keywords}.R"
	script:
		"scripts/get_samples.R"

# ------------------------------------------------------------------------------
# Extract and process expression data for specific keywords from main h5 file
# ------------------------------------------------------------------------------
rule extract_data:
	input:
		h5="results/human_matrix_download.h5",
		samples="results/sample_definitions/{keywords}.R"
	output:
		expr=config["data_dir"] + "{keywords}/expression_matrix_norm_sva.tsv",
		raw=config["data_dir"] + "{keywords}/expression_matrix_raw.tsv",
		design=config["data_dir"] + "{keywords}/design.tsv"
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
		expr=config["data_dir"] + "{keywords}/expression_matrix_norm_sva.tsv",
		raw=config["data_dir"] + "{keywords}/expression_matrix_raw.tsv",
		design=config["data_dir"] + "{keywords}/design.tsv"
	output:
		config["data_dir"] + "{keywords}/summary.html"
	params:
		pdf_out = config["data_dir"] + "{keywords}/tsne_plots.pdf"
	script:
		"scripts/explore_data.Rmd"

# ------------------------------------------------------------------------------
# Target rule to process some interesting tissues
# ------------------------------------------------------------------------------
rule explore_thyroid_blood_muscle:
	input:
		config["data_dir"] + "Thyroid/summary.html",
		config["data_dir"] + "Whole_Blood/summary.html",
		config["data_dir"] + "Skeletal_Muscle/summary.html"
	output:
		"results/summaries.zip"
	shell:
		"zip {output[0]} {input}"
