# Configuration

You need one configuration file and one annotation file to run the complete workflow. Additionally, you can provide a feature annotation file in the project configuration (e.g., for plotting gene symbols instead of ensembl terms). If in doubt read the comments in the config and/or try the default values.

- project configuration (`config/config.yaml`): configures the analyses to be performed and is different for every project/dataset.
- annotation (annotation): CSV file with one row per analysis and consisting of 10 mandatory columns
    -  name: name of the dataset/analysis (tip: keep it short, but descriptive, distinctive and unique; `OvA` is a forbidden keyword)
    -  data: Absolute path to the input data as CSV file as feature by sample table (eg RNA count matrix) that has already been quality controlled (eg bad samples removed) and filtered for relevant features (eg only expressed genes). The first column has to contain the features and the first row the sample-names.
    -  metadata: Absolute path to the metadata as CSV file for the required analysis (ie every variable in the formula and optional blocking variable needs to have a corresponding column). The first column has to be the sample name. The metadata file has to be R compatible (eg column names should not start with a number or contain colons).
    -  formula: A string that will be converted to a formula in R (eg ~ treatment + batch).
    -  block_var: Flag to indicate which variable (present in the metadata as column) should be used for the blocking feature (see README > Features) or if it should be skipped (0).
    -  comparisons: Variable names contained in the formula (and metadata) which coefficient's you are interested in, separated by `|` (e.g., `treatment|batch`). Results of all derived groups (e.g., `treatmentLPS`) containing one of the comparisons will be returned in the results. Empty comparisons field is allowed if one is only interested in the one-vs-all contrasts provided by the OvA analysis (see OvA field) and results in empty result files and plots.
    -  calcNormFactors_method: Flag to indicate if [edgeR:calcNormFactors](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/calcNormFactors) function should be used specifing the parameter "method" (eg none or TMM) or should be skipped because the input data is already log-normalized (0).
    -  voom: Flag to indicate if voom function should be used (1) or not (0). Note: Should be 0 if data is already log-normalized and/or limma-trend will be used.
    -  eBayes: Flag to indicate if eBayes function should be used (1) or not (0). Note: Skipping eBayes (0) will lead to the use of ordinary t-statistic with topTable and is [not recommended by the limma author Gordon Smyth](https://support.bioconductor.org/p/35174/), the B-statistics (log-odds) are still determined using eBayes, assuming they will not be used downstream. Make sure you know what you are doing.
    -  limma_trend: Flag to indicate if limma-trend should be used (1) (ie sets [limma::eBayes](https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/ebayes) parameter trend=TRUE), or not (0). Note: Make sure to activate the required eBayes function (=1) and deactivate voom (=0) if you use limma-trend. Using both, voom and limma-trend makes no sense, but is not prohibited by the workflow.
    -  OvA (optional): Variable names contained in the formula (and metadata) for which a one-vs-all analysis using contrasts should be performed, separated by `|` (e.g., `cell_type|batch`). OvA for interaction terms and continuous variables is not supported and causes an error. Formulas that contain interaction terms are allowed, but should be used with care and only if questions cannot be addressed with simpler models, since they make both main and interaction effects' OvA results change with reference levels. Due to numerical approximations used by limma [dicussed here](https://support.bioconductor.org/p/112788/), p-values change depending on whether a term is modeled as means model (`~ 0 + term1_as_means_model + term2_as_mean_ref_model + term3_as_mean_ref_model`) or as mean-reference model (`~ term1_as_mean_ref_model + term2_as_mean_ref_model + term3_as_mean_ref_model`) - the logFCs are stable. Keep this in mind for the interpretation, i.e., be careful about genes that have p-values close to the significance threshold and low logFCs. Otherwise, the order of terms does not matter.

Set workflow-specific `resources` or command line arguments (CLI) in the workflow profile `workflow/profiles/default.config.yaml`, which supersedes global Snakemake profiles.

## Common configuration scenarios
- standard **limma-voom** workflow with raw counts as input data ([see "Differential expression: voom" in the limma userguide](http://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf))
    - calcNormFactors_method: none (or other normalization method e.g., TMM, to be considered by voom)
    - voom: 1
    - eBayes: 1
    - limma_trend: 0
- standard **limma-trend** workflow with raw counts as input data ([see "Differential expression: limma-trend" in the limma userguide](http://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf))
    - calcNormFactors_method: none (or other normalization method e.g., TMM)
    - voom: 0
    - eBayes: 1
    - limma_trend: 1
- limma-trend workflow for **log-normalized input data**
    - calcNormFactors_method: 0 (thereby skipping any normalization or voom related actions)
    - voom: 0
    - eBayes: 1
    - limma_trend: 1
- for more in-depth understanding check out the commented code: [limma.R](../workflow/scripts/limma.R)

