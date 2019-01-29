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
# Download and process expression data for specific keywords
# ------------------------------------------------------------------------------
rule download_samples:
	input:
		h5="results/human_matrix_download.h5"
	output:
		expr="results/downloads/{keywords}/expression_normalized.tsv",
		design="results/downloads/{keywords}/design.tsv"
	log:
		"logs/download_samples_{keywords}.log"
	conda:
		"envs/r_env.yaml"
	params:
		norm_method = "sva" # must be one of 'sva', 'peer' or 'quantile'
	script:
		"scripts/download.R"

# ------------------------------------------------------------------------------
# Explore the data and create a nice summary
# ------------------------------------------------------------------------------
rule explore_data:
	input:
		expr="results/downloads/{keywords}/expression_normalized.tsv",
		design="results/downloads/{keywords}/design.tsv"
	output:
		"results/downloads/{keywords}/summary.html"
	conda:
		"envs/r_env.yaml"
	wildcard_constraints:
		keywords = "^(?!all_samples)$"
	script:
		"scripts/explore_data.Rmd"
