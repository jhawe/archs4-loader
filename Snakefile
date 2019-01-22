rule create_sample_sheet:
	output:
		"results/design_all_samples.tsv",
		"results/human_matrix_download.h5"
	script:
		"scripts/create_sample_sheet.R"

rule get_samples:
	input:
		design="results/design_all_samples.tsv",
	output:
		"results/sample_definitions/{keywords}.R"
	script:
		"scripts/get_samples.R"

rule download_samples:
	input:
		h5="results/human_matrix_download.h5",
		samples="results/sample_definitions/{keywords}.R"
	output:
		expr="results/downloads/{keywords}/expression_matrix_norm_sva.tsv",
		plot="results/downloads/{keywords}/expression_matrix_norm_sva.pdf",
		design="results/downloads/{keywords}/design.tsv"
	log:
		"logs/download_samples_{keywords}.log"
	script:
		"scripts/download.R"

rule explore_data:
	input:
		expr="results/downloads/{keywords}/expression_matrix_norm_sva.tsv"
	output:
		"results/downloads/{keywords}/exploration.pdf"
	script:
		"scripts/explore_data.Rmd"
