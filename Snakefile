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
# Download and process expression data for specific keywords
# ------------------------------------------------------------------------------
rule download_samples:
	input:
		h5="results/human_matrix_download.h5",
		samples="results/sample_definitions/{keywords}.R"
	output:
		expr="results/downloads/{keywords}/expression_matrix_norm_sva.tsv",
#		plot="results/downloads/{keywords}/expression_matrix_norm_sva.pdf",
		design="results/downloads/{keywords}/design.tsv"
	log:
		"logs/download_samples_{keywords}.log"
	script:
		"scripts/download.R"

# ------------------------------------------------------------------------------
# Explore the data and create a nice summary
# ------------------------------------------------------------------------------
rule explore_data:
	input:
		expr="results/downloads/{keywords}/expression_matrix_norm_sva.tsv",
		design="results/downloads/{keywords}/design.tsv"
	output:
		"results/downloads/{keywords}/summary.html"
	script:
		"scripts/explore_data.Rmd"
