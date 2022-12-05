---
title: "Run differential gene expression and gene-sets analysis"
author: "Francesc Castro-Giner"
date: "2022-02-23"
output: workflowr::wflow_html
editor_options:
  chunk_output_type: console
params:
  date: 'October 13, 2022'
  sce_dir: ./data/sce
  min_counts: 5
  min_present_prop: 0.50
---

<p>
<button type="button" class="btn btn-default btn-workflowr btn-workflowr-report"
  data-toggle="collapse" data-target="#workflowr-report">
  <span class="glyphicon glyphicon-list" aria-hidden="true"></span>
  workflowr
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
</button>
</p>

<div id="workflowr-report" class="collapse">
<ul class="nav nav-tabs">
  <li class="active"><a data-toggle="tab" href="#summary">Summary</a></li>
  <li><a data-toggle="tab" href="#checks">
  Checks <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  </a></li>
  <li><a data-toggle="tab" href="#versions">Past versions</a></li>
</ul>

<div class="tab-content">
<div id="summary" class="tab-pane fade in active">
  <p><strong>Last updated:</strong> 2022-10-13</p>
  <p><strong>Checks:</strong>
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  7
  <span class="glyphicon glyphicon-exclamation-sign text-danger" aria-hidden="true"></span>
  0
  </p>
  <p><strong>Knit directory:</strong>
  <code>diamantopoulou-ctc-dynamics/</code>
  <span class="glyphicon glyphicon-question-sign" aria-hidden="true"
  title="This is the local directory in which the code in this file was executed.">
  </span>
  </p>
  <p>
  This reproducible <a href="https://rmarkdown.rstudio.com">R Markdown</a>
  analysis was created with <a
  href="https://github.com/workflowr/workflowr">workflowr</a> (version
  1.7.0). The <em>Checks</em> tab describes the
  reproducibility checks that were applied when the results were created.
  The <em>Past versions</em> tab lists the development history.
  </p>
<hr>
</div>
<div id="checks" class="tab-pane fade">
  <div class="panel-group" id="workflowr-checks">
  <div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongRMarkdownfilestronguptodate">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>R Markdown file:</strong> up-to-date
</a>
</p>
</div>
<div id="strongRMarkdownfilestronguptodate" class="panel-collapse collapse">
<div class="panel-body">
  
Great! Since the R Markdown file has been committed to the Git repository, you
know the exact version of the code that produced these results.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongEnvironmentstrongempty">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Environment:</strong> empty
</a>
</p>
</div>
<div id="strongEnvironmentstrongempty" class="panel-collapse collapse">
<div class="panel-body">
  
Great job! The global environment was empty. Objects defined in the global
environment can affect the analysis in your R Markdown file in unknown ways.
For reproduciblity it's best to always run the code in an empty environment.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongSeedstrongcodesetseed20220425code">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Seed:</strong> <code>set.seed(20220425)</code>
</a>
</p>
</div>
<div id="strongSeedstrongcodesetseed20220425code" class="panel-collapse collapse">
<div class="panel-body">
  
The command <code>set.seed(20220425)</code> was run prior to running the code in the R Markdown file.
Setting a seed ensures that any results that rely on randomness, e.g.
subsampling or permutations, are reproducible.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongSessioninformationstrongrecorded">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Session information:</strong> recorded
</a>
</p>
</div>
<div id="strongSessioninformationstrongrecorded" class="panel-collapse collapse">
<div class="panel-body">
  
Great job! Recording the operating system, R version, and package versions is
critical for reproducibility.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongCachestrongnone">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Cache:</strong> none
</a>
</p>
</div>
<div id="strongCachestrongnone" class="panel-collapse collapse">
<div class="panel-body">
  
Nice! There were no cached chunks for this analysis, so you can be confident
that you successfully produced the results during this run.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongFilepathsstrongrelative">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>File paths:</strong> relative
</a>
</p>
</div>
<div id="strongFilepathsstrongrelative" class="panel-collapse collapse">
<div class="panel-body">
  
Great job! Using relative paths to the files within your workflowr project
makes it easier to run your code on other machines.

</div>
</div>
</div>
<div class="panel panel-default">
<div class="panel-heading">
<p class="panel-title">
<a data-toggle="collapse" data-parent="#workflowr-checks" href="#strongRepositoryversionstrongahrefhttpsgithubcomVUQCBdiamantopoulouctcdynamicstree035a69de313bceade92e9e94e6d61311e1fe4288targetblank035a69da">
  <span class="glyphicon glyphicon-ok text-success" aria-hidden="true"></span>
  <strong>Repository version:</strong> <a href="https://github.com/VU-QCB/diamantopoulou-ctc-dynamics/tree/035a69de313bceade92e9e94e6d61311e1fe4288" target="_blank">035a69d</a>
</a>
</p>
</div>
<div id="strongRepositoryversionstrongahrefhttpsgithubcomVUQCBdiamantopoulouctcdynamicstree035a69de313bceade92e9e94e6d61311e1fe4288targetblank035a69da" class="panel-collapse collapse">
<div class="panel-body">
  
<p>
Great! You are using Git for version control. Tracking code development and
connecting the code version to the results is critical for reproducibility.
</p>

<p>
The results in this page were generated with repository version <a href="https://github.com/VU-QCB/diamantopoulou-ctc-dynamics/tree/035a69de313bceade92e9e94e6d61311e1fe4288" target="_blank">035a69d</a>.
See the <em>Past versions</em> tab to see a history of the changes made to the
R Markdown and HTML files.
</p>

<p>
Note that you need to be careful to ensure that all relevant files for the
analysis have been committed to Git prior to generating the results (you can
use <code>wflow_publish</code> or <code>wflow_git_commit</code>). workflowr only
checks the R Markdown file, but you know if there are other scripts or data
files that it depends on. Below is the status of the Git repository when the
results were generated:
</p>

<pre><code>
Ignored files:
	Ignored:    .DS_Store
	Ignored:    .Rproj.user/
	Ignored:    code/.DS_Store
	Ignored:    configuration/.DS_Store
	Ignored:    data/.DS_Store
	Ignored:    data/differential_expression/
	Ignored:    data/patients/
	Ignored:    data/resources/
	Ignored:    data/sce/

Unstaged changes:
	Modified:   code/R-functions/install_reqd_pkgs.r

</code></pre>

<p>
Note that any generated files, e.g. HTML, png, CSS, etc., are not included in
this status report because it is ok for generated content to have uncommitted
changes.
</p>

</div>
</div>
</div>
</div>
<hr>
</div>
<div id="versions" class="tab-pane fade">
  
<p>
These are the previous versions of the repository in which changes were made
to the R Markdown (<code>analysis/0_differential_expression_gsea_gsva.Rmd</code>) and HTML (<code>docs/0_differential_expression_gsea_gsva.html</code>)
files. If you've configured a remote Git repository (see
<code>?wflow_git_remote</code>), click on the hyperlinks in the table below to
view the files as they were in that past version.
</p>
<div class="table-responsive">
<table class="table table-condensed table-hover">
<thead>
<tr>
<th>File</th>
<th>Version</th>
<th>Author</th>
<th>Date</th>
<th>Message</th>
</tr>
</thead>
<tbody>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/VU-QCB/diamantopoulou-ctc-dynamics/41453ae72586bc287e8d7d9e634aa81822737479/docs/0_differential_expression_gsea_gsva.html" target="_blank">41453ae</a></td>
<td>fcg-bio</td>
<td>2022-05-12</td>
<td>Build site.</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/VU-QCB/diamantopoulou-ctc-dynamics/1fc87b5937dd4e880b8ee9f151f8f8ab6a5c0b3a/docs/0_differential_expression_gsea_gsva.html" target="_blank">1fc87b5</a></td>
<td>fcg-bio</td>
<td>2022-05-10</td>
<td>Build site.</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/VU-QCB/diamantopoulou-ctc-dynamics/8fb5513f546f6646b28b879c3611e491e358b4cd/docs/0_differential_expression_gsea_gsva.html" target="_blank">8fb5513</a></td>
<td>fcg-bio</td>
<td>2022-05-10</td>
<td>Build site.</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/VU-QCB/diamantopoulou-ctc-dynamics/74b189131a39d49fd9e79eaee6c8e76575069de6/docs/0_differential_expression_gsea_gsva.html" target="_blank">74b1891</a></td>
<td>fcg-bio</td>
<td>2022-04-26</td>
<td>Build site.</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/VU-QCB/diamantopoulou-ctc-dynamics/c0865c6165c3d4748ca1d329af7ee9a649bc865c/docs/0_differential_expression_gsea_gsva.html" target="_blank">c0865c6</a></td>
<td>fcg-bio</td>
<td>2022-04-26</td>
<td>Build site.</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/VU-QCB/diamantopoulou-ctc-dynamics/a1365909e9e1773494c33418c8baf5adb462c53d/docs/0_differential_expression_gsea_gsva.html" target="_blank">a136590</a></td>
<td>fcg-bio</td>
<td>2022-04-26</td>
<td>Build site.</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/VU-QCB/diamantopoulou-ctc-dynamics/bfb622b5aac934265e71dc66dd5a0985a4fca3d1/docs/0_differential_expression_gsea_gsva.html" target="_blank">bfb622b</a></td>
<td>fcg-bio</td>
<td>2022-04-26</td>
<td>Build site.</td>
</tr>
<tr>
<td>html</td>
<td><a href="https://rawcdn.githack.com/VU-QCB/diamantopoulou-ctc-dynamics/1006c84c731de8ffc432a602cb2472a755530e58/docs/0_differential_expression_gsea_gsva.html" target="_blank">1006c84</a></td>
<td>fcg-bio</td>
<td>2022-04-25</td>
<td>Build site.</td>
</tr>
<tr>
<td>Rmd</td>
<td><a href="https://github.com/VU-QCB/diamantopoulou-ctc-dynamics/blob/0ded9f5eb929c7853186f31da3ed481e57e151e9/analysis/0_differential_expression_gsea_gsva.Rmd" target="_blank">0ded9f5</a></td>
<td>fcg-bio</td>
<td>2022-04-25</td>
<td>added final code</td>
</tr>
</tbody>
</table>
</div>

<hr>
</div>
</div>
</div>






## Load libraries, additional functions and data

Setup environment

```r
knitr::opts_chunk$set(results='asis', echo=TRUE, message=FALSE, warning=FALSE, error=FALSE, fig.align = 'center', fig.width = 3.5, fig.asp = 0.618, dpi = 600, dev = c("png", "pdf"), fig.showtext = TRUE)

options(stringsAsFactors = FALSE)
```

Load packages

```r
library(tidyverse)
library(scater)
library(scran)
library(edgeR)
library(clusterProfiler)
library(GSVA)
library(foreach)
```

Load shared variables

```r
source("./configuration/rmarkdown/shared_variables.R")
```

Load custom functions

```r
source('./code/R-functions/dge_wrappers.r')
source('./code/R-functions/gse_omnibus.r')
source('./code/R-functions/gse_report.r')
clean_msigdb_names <- function(x) x %>% gsub('REACTOME_', '', .) %>% gsub('WP_', '', .) %>% gsub('BIOCARTA_', '', .) %>% gsub('KEGG_', '', .) %>% gsub('PID_', '', .) %>% gsub('GOBP_', '', .) %>% gsub('_', ' ', .)
```

Load MSigDB gene sets

```r
gmt_files_symbols <- list(
  msigdb.c2.cp = './data/resources/MSigDB/v7.4/c2.cp.v7.4.symbols.gmt',
  msigdb.c5.bp = './data/resources/MSigDB/v7.4/c5.go.bp.v7.4.symbols.gmt'
)

gmt_files_entrez <- list(
  msigdb.c2.cp = './data/resources/MSigDB/v7.4/c2.cp.v7.4.entrez.gmt',
  msigdb.c5.bp = './data/resources/MSigDB/v7.4/c5.go.bp.v7.4.entrez.gmt'
)

# combine MSigDB.C2.CP and GO:BP
new_file <- gsub('c2.cp', 'c2.cp.c5.bp', gmt_files_symbols$msigdb.c2.cp)
cat_cmd <- paste('cat', gmt_files_symbols$msigdb.c5.bp,  gmt_files_symbols$msigdb.c2.cp, '>',new_file)
system(cat_cmd)
gmt_files_symbols$msigdb.c2.cp.c5.bp <- new_file

gmt_sets <- lapply(gmt_files_symbols, function(x) read.gmt(x) %>% collect %>% .[['term']] %>% levels)
```

## NSG-CDX-BR16 : all samples
Configuration

```r
use_sce <- readRDS(file = file.path(params$sce_dir, 'sce_br16.rds'))
output_dir <- './data/differential_expression/br16'
if(!file.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
```

Run DGE analysis

```r
dge <- edgeR_dge(
  use_sce,
  # Desing configuration for differential expression
  group_var =  'timepoint',
  group_sample = 'resting',
  group_ref = 'active',
  numeric_covar = NULL,
  batch_vars = NULL,
  design_formula = "~ 0 + timepoint",
  coef = 'last',
  # Conversion from SingleCellExperiment to DGEList
  spike_normalization = FALSE,
  assay_to_DGEList = 'counts',
  assay_to_row_filter = "counts",
  use_colData = NULL,
  use_rowData = NULL,
  # Feature filtering parameters
  use_filterByExpr = TRUE,
  min_counts = params$min_counts,
  min_present_prop = params$min_present_prop,
  # EdgeR workflow configuration
  run_calcNormFactors = 'TMM',
  estimateDisp_robust = FALSE,
  estimateDisp_trend.method = "locfit",
  glmQLFit_robust = TRUE,
  glm_approach = "QLF",
  # Output configuration
  adjust_method = 'BH',
  assays_from_SingleCellExperiment = NULL
  )

# Add gene description
httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl <-  biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_desc <- biomaRt::getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = dge$results$gene_name, mart =ensembl) %>% 
  dplyr::rename('gene_name' = 'external_gene_name')
use_res <- dge$results %>%  left_join(., gene_desc)
dge$results <- use_res %>% 
  filter(!duplicated(feature)) %>% 
  mutate(rownames = feature) %>% 
  column_to_rownames('rownames')

detach("package:biomaRt", unload=TRUE)

saveRDS(dge, file = file.path(output_dir, 'dge_edgeR_QLF_robust.rds'))
```

Run GSEA

```r
dge <- readRDS(file.path(output_dir, 'dge_edgeR_QLF_robust.rds'))
res_gse <- gse_omnibus(
    feature_names = dge$results$gene_name,
    p = dge$results$FDR,
    fc = dge$results$logFC,
    gmt_files = gmt_files_symbols, 

    save_intermediates = file.path(output_dir, 'gse_omnibus'),
    
    run_all_ora = FALSE,
    run_all_gsea = FALSE,
    run_GSEA = TRUE,
    run_gseGO = FALSE,

    args_gse = list(minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1),

    )
saveRDS(res_gse, file = file.path(output_dir, 'gse_gsea.rds'))
```

Clean data

```r
rm(use_sce)
rm(dge)
rm(res_gse)
```

## NSG-CDX-BR16 : CTC-Cluster and CTC-WBC
Configuration

```r
use_sce <- use_sce[,use_sce$sample_type_g == 'ctc_cluster']
output_dir <- './data/differential_expression/br16-ctc_cluster_and_wbc'
if(!file.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
```

Run DGE analysis

```r
dge <- edgeR_dge(
  use_sce,
  # Desing configuration for differential expression
  group_var =  'timepoint',
  group_sample = 'resting',
  group_ref = 'active',
  numeric_covar = NULL,
  batch_vars = NULL,
  design_formula = "~ 0 + timepoint",
  coef = 'last',
  # Conversion from SingleCellExperiment to DGEList
  spike_normalization = FALSE,
  assay_to_DGEList = 'counts',
  assay_to_row_filter = "counts",
  use_colData = NULL,
  use_rowData = NULL,
  # Feature filtering parameters
  use_filterByExpr = TRUE,
  min_counts = params$min_counts,
  min_present_prop = params$min_present_prop,
  # EdgeR workflow configuration
  run_calcNormFactors = 'TMM',
  estimateDisp_robust = FALSE,
  estimateDisp_trend.method = "locfit",
  glmQLFit_robust = TRUE,
  glm_approach = "QLF",
  # Output configuration
  adjust_method = 'BH',
  assays_from_SingleCellExperiment = NULL
  )

# Add gene description
httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl <-  biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_desc <- biomaRt::getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = dge$results$gene_name, mart =ensembl) %>% 
  dplyr::rename('gene_name' = 'external_gene_name')
use_res <- dge$results %>%  left_join(., gene_desc)
dge$results <- use_res %>% 
  filter(!duplicated(feature)) %>% 
  mutate(rownames = feature) %>% 
  column_to_rownames('rownames')

detach("package:biomaRt", unload=TRUE)

saveRDS(dge, file = file.path(output_dir, 'dge_edgeR_QLF_robust.rds'))
```

Run GSEA

```r
dge <- readRDS(file.path(output_dir, 'dge_edgeR_QLF_robust.rds'))
res_gse <- gse_omnibus(
    feature_names = dge$results$gene_name,
    p = dge$results$FDR,
    fc = dge$results$logFC,
    gmt_files = gmt_files_symbols, 

    save_intermediates = file.path(output_dir, 'gse_omnibus'),
    
    run_all_ora = FALSE,
    run_all_gsea = FALSE,
    run_GSEA = TRUE,
    run_gseGO = FALSE,

    args_gse = list(minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1),

    )
saveRDS(res_gse, file = file.path(output_dir, 'gse_gsea.rds'))
```

Clean data

```r
rm(use_sce)
rm(dge)
rm(res_gse)
```


## NSG-CDX-BR16 : CTC-Single
Configuration

```r
use_sce <- use_sce[,use_sce$sample_type_g == 'ctc_single']
output_dir <- './data/differential_expression/br16-ctc_single'
if(!file.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
```

Run DGE analysis

```r
dge <- edgeR_dge(
  use_sce,
  # Desing configuration for differential expression
  group_var =  'timepoint',
  group_sample = 'resting',
  group_ref = 'active',
  numeric_covar = NULL,
  batch_vars = NULL,
  design_formula = "~ 0 + timepoint",
  coef = 'last',
  # Conversion from SingleCellExperiment to DGEList
  spike_normalization = FALSE,
  assay_to_DGEList = 'counts',
  assay_to_row_filter = "counts",
  use_colData = NULL,
  use_rowData = NULL,
  # Feature filtering parameters
  use_filterByExpr = TRUE,
  min_counts = params$min_counts,
  min_present_prop = params$min_present_prop,
  # EdgeR workflow configuration
  run_calcNormFactors = 'TMM',
  estimateDisp_robust = FALSE,
  estimateDisp_trend.method = "locfit",
  glmQLFit_robust = TRUE,
  glm_approach = "QLF",
  # Output configuration
  adjust_method = 'BH',
  assays_from_SingleCellExperiment = NULL
  )

# Add gene description
httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl <-  biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_desc <- biomaRt::getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = dge$results$gene_name, mart =ensembl) %>% 
  dplyr::rename('gene_name' = 'external_gene_name')
use_res <- dge$results %>%  left_join(., gene_desc)
dge$results <- use_res %>% 
  filter(!duplicated(feature)) %>% 
  mutate(rownames = feature) %>% 
  column_to_rownames('rownames')

detach("package:biomaRt", unload=TRUE)

saveRDS(dge, file = file.path(output_dir, 'dge_edgeR_QLF_robust.rds'))
```

Run GSEA

```r
dge <- readRDS(file.path(output_dir, 'dge_edgeR_QLF_robust.rds'))
res_gse <- gse_omnibus(
    feature_names = dge$results$gene_name,
    p = dge$results$FDR,
    fc = dge$results$logFC,
    gmt_files = gmt_files_symbols, 

    save_intermediates = file.path(output_dir, 'gse_omnibus'),
    
    run_all_ora = FALSE,
    run_all_gsea = FALSE,
    run_GSEA = TRUE,
    run_gseGO = FALSE,

    args_gse = list(minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1),

    )
saveRDS(res_gse, file = file.path(output_dir, 'gse_gsea.rds'))
```

Clean data

```r
rm(use_sce)
rm(dge)
rm(res_gse)
```




## NSG-LM2
Configuration

```r
use_sce <- readRDS(file = file.path(params$sce_dir, 'sce_lm2.rds'))
output_dir <- './data/differential_expression/lm2'
if(!file.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
```

Run DGE analysis

```r
dge <- edgeR_dge(
  use_sce,
  # Desing configuration for differential expression
  group_var =  'timepoint',
  group_sample = 'resting',
  group_ref = 'active',
  numeric_covar = NULL,
  batch_vars = NULL,
  design_formula = "~ 0 + timepoint",
  coef = 'last',
  # Conversion from SingleCellExperiment to DGEList
  spike_normalization = FALSE,
  assay_to_DGEList = 'counts',
  assay_to_row_filter = "counts",
  use_colData = NULL,
  use_rowData = NULL,
  # Feature filtering parameters
  use_filterByExpr = TRUE,
  min_counts = params$min_counts,
  min_present_prop = params$min_present_prop,
  # EdgeR workflow configuration
  run_calcNormFactors = 'TMM',
  estimateDisp_robust = FALSE,
  estimateDisp_trend.method = "locfit",
  glmQLFit_robust = TRUE,
  glm_approach = "QLF",
  # Output configuration
  adjust_method = 'BH',
  assays_from_SingleCellExperiment = NULL
  )

# Add gene description
httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl <-  biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_desc <- biomaRt::getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = dge$results$gene_name, mart =ensembl) %>% 
  dplyr::rename('gene_name' = 'external_gene_name')
use_res <- dge$results %>%  left_join(., gene_desc)
dge$results <- use_res %>% 
  filter(!duplicated(feature)) %>% 
  mutate(rownames = feature) %>% 
  column_to_rownames('rownames')

detach("package:biomaRt", unload=TRUE)

saveRDS(dge, file = file.path(output_dir, 'dge_edgeR_QLF_robust.rds'))
```

Clean data

```r
rm(use_sce)
rm(dge)
```


## Patient
Configuration

```r
use_sce <- readRDS(file = file.path(params$sce_dir, 'sce_patient.rds'))
output_dir <- './data/differential_expression/patient'
if(!file.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
```

Run DGE analysis

```r
dge <- edgeR_dge(
  use_sce,
  # Desing configuration for differential expression
  group_var =  'timepoint',
  group_sample = 'resting',
  group_ref = 'active',
  numeric_covar = NULL,
  batch_vars = NULL,
  design_formula = "~ 0 + timepoint",
  coef = 'last',
  # Conversion from SingleCellExperiment to DGEList
  spike_normalization = FALSE,
  assay_to_DGEList = 'counts',
  assay_to_row_filter = "counts",
  use_colData = NULL,
  use_rowData = NULL,
  # Feature filtering parameters
  use_filterByExpr = TRUE,
  min_counts = params$min_counts,
  min_present_prop = params$min_present_prop,
  # EdgeR workflow configuration
  run_calcNormFactors = 'TMM',
  estimateDisp_robust = FALSE,
  estimateDisp_trend.method = "locfit",
  glmQLFit_robust = TRUE,
  glm_approach = "QLF",
  # Output configuration
  adjust_method = 'BH',
  assays_from_SingleCellExperiment = NULL
  )

# Add gene description
httr::set_config(httr::config(ssl_verifypeer = FALSE))
ensembl <-  biomaRt::useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
gene_desc <- biomaRt::getBM(attributes=c('external_gene_name','description'), filters = 'external_gene_name', values = dge$results$gene_name, mart =ensembl) %>% 
  dplyr::rename('gene_name' = 'external_gene_name')
use_res <- dge$results %>%  left_join(., gene_desc)
dge$results <- use_res %>% 
  filter(!duplicated(feature)) %>% 
  mutate(rownames = feature) %>% 
  column_to_rownames('rownames')

detach("package:biomaRt", unload=TRUE)

saveRDS(dge, file = file.path(output_dir, 'dge_edgeR_QLF_robust.rds'))
```

Clean data

```r
rm(use_sce)
rm(dge)
```

## LM2 time kinetics
Configuration

```r
use_sce <- readRDS(file = file.path(params$sce_dir, 'sce_lm2_tk.rds'))
output_dir <- './data/differential_expression/lm2_tk'
if(!file.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)
```

Run [GSVA](https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html) run with gene-set size between 5 and 700. Original GSEA analysis was performed with 10-500, but with this new treshold we make sure that all the gene sets from BR16 results are included in the analysis, as the effective gene set (expressed genes) might be different in GSVA analysis.

For this analysis we remove samples from timepoint ZT0 (06:00). It only contains one replicate and can bias results. The timepoint will be added for visualization.

```r
use_sce <- use_sce[,!use_sce$timepoint %in% c('0600')]
rownames(use_sce) <- rowData(use_sce)$gene_name
use_gmt_file <- "./data/resources/MSigDB/v7.4/c2.cp.c5.bp.v7.4.symbols.gmt"
gset <- GSEABase::getGmt(use_gmt_file)
gset_db <- foreach(x = gset, .combine = rbind) %do% {c(term_size = length(x@geneIds))} %>% data.frame()
gset_db$term_name <- names(gset)

gsva_res <- gsva(assay(use_sce, 'logcpm'), 
                   method = 'gsva',
                   gset.idx.list = gset, 
                   min.sz = 5, 
                   max.sz = 700, 
                   kcdf = "Gaussian",
                   mx.diff = TRUE, 
                   verbose = FALSE)

saveRDS(gsva_res, file = file.path(output_dir, 'gsva_c2.cp.c5.bp.rds'))
```




<br>
<p>
<button type="button" class="btn btn-default btn-workflowr btn-workflowr-sessioninfo"
  data-toggle="collapse" data-target="#workflowr-sessioninfo"
  style = "display: block;">
  <span class="glyphicon glyphicon-wrench" aria-hidden="true"></span>
  Session information
</button>
</p>

<div id="workflowr-sessioninfo" class="collapse">

```r
sessionInfo()
```

R version 4.2.1 (2022-06-23)
Platform: aarch64-apple-darwin20 (64-bit)
Running under: macOS Monterey 12.6

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods  
[8] base     

other attached packages:
 [1] foreach_1.5.2               GSVA_1.44.5                
 [3] clusterProfiler_4.4.4       edgeR_3.38.4               
 [5] limma_3.52.4                scran_1.24.1               
 [7] scater_1.24.0               scuttle_1.6.3              
 [9] SingleCellExperiment_1.18.1 SummarizedExperiment_1.26.1
[11] Biobase_2.56.0              GenomicRanges_1.48.0       
[13] GenomeInfoDb_1.32.4         IRanges_2.30.1             
[15] S4Vectors_0.34.0            BiocGenerics_0.42.0        
[17] MatrixGenerics_1.8.1        matrixStats_0.62.0         
[19] forcats_0.5.2               stringr_1.4.1              
[21] dplyr_1.0.10                purrr_0.3.5                
[23] readr_2.1.3                 tidyr_1.2.1                
[25] tibble_3.1.8                ggplot2_3.3.6              
[27] tidyverse_1.3.2            

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                tidyselect_1.2.0         
  [3] RSQLite_2.2.18            AnnotationDbi_1.58.0     
  [5] grid_4.2.1                BiocParallel_1.30.4      
  [7] scatterpie_0.1.8          munsell_0.5.0            
  [9] ScaledMatrix_1.4.1        codetools_0.2-18         
 [11] statmod_1.4.37            withr_2.5.0              
 [13] colorspace_2.0-3          GOSemSim_2.22.0          
 [15] knitr_1.40                rstudioapi_0.14          
 [17] DOSE_3.22.1               git2r_0.30.1             
 [19] GenomeInfoDbData_1.2.8    polyclip_1.10-0          
 [21] bit64_4.0.5               farver_2.1.1             
 [23] rhdf5_2.40.0              rprojroot_2.0.3          
 [25] downloader_0.4            treeio_1.20.2            
 [27] vctrs_0.4.2               generics_0.1.3           
 [29] xfun_0.33                 R6_2.5.1                 
 [31] ggbeeswarm_0.6.0          graphlayouts_0.8.2       
 [33] rsvd_1.0.5                locfit_1.5-9.6           
 [35] rhdf5filters_1.8.0        bitops_1.0-7             
 [37] cachem_1.0.6              fgsea_1.22.0             
 [39] gridGraphics_0.5-1        DelayedArray_0.22.0      
 [41] assertthat_0.2.1          showtext_0.9-5           
 [43] promises_1.2.0.1          scales_1.2.1             
 [45] ggraph_2.1.0              enrichplot_1.16.2        
 [47] googlesheets4_1.0.1       beeswarm_0.4.0           
 [49] gtable_0.3.1              beachmat_2.12.0          
 [51] tidygraph_1.2.2           workflowr_1.7.0          
 [53] rlang_1.0.6               splines_4.2.1            
 [55] lazyeval_0.2.2            gargle_1.2.1             
 [57] broom_1.0.1               yaml_2.3.5               
 [59] reshape2_1.4.4            modelr_0.1.9             
 [61] backports_1.4.1           httpuv_1.6.6             
 [63] qvalue_2.28.0             tools_4.2.1              
 [65] ggplotify_0.1.0           ellipsis_0.3.2           
 [67] jquerylib_0.1.4           RColorBrewer_1.1-3       
 [69] Rcpp_1.0.9                plyr_1.8.7               
 [71] sparseMatrixStats_1.8.0   zlibbioc_1.42.0          
 [73] RCurl_1.98-1.9            viridis_0.6.2            
 [75] haven_2.5.1               ggrepel_0.9.1            
 [77] cluster_2.1.4             fs_1.5.2                 
 [79] magrittr_2.0.3            data.table_1.14.2        
 [81] DO.db_2.9                 reprex_2.0.2             
 [83] googledrive_2.0.0         whisker_0.4              
 [85] xtable_1.8-4              hms_1.1.2                
 [87] patchwork_1.1.2           evaluate_0.17            
 [89] XML_3.99-0.11             readxl_1.4.1             
 [91] gridExtra_2.3             compiler_4.2.1           
 [93] shadowtext_0.1.2          crayon_1.5.2             
 [95] htmltools_0.5.3           ggfun_0.0.7              
 [97] later_1.3.0               tzdb_0.3.0               
 [99] aplot_0.1.8               lubridate_1.8.0          
[101] DBI_1.1.3                 tweenr_2.0.2             
[103] dbplyr_2.2.1              MASS_7.3-58.1            
[105] Matrix_1.5-1              cli_3.4.1                
[107] parallel_4.2.1            metapod_1.4.0            
[109] igraph_1.3.5              pkgconfig_2.0.3          
[111] xml2_1.3.3                annotate_1.74.0          
[113] ggtree_3.4.4              vipor_0.4.5              
[115] bslib_0.4.0               dqrng_0.3.0              
[117] XVector_0.36.0            rvest_1.0.3              
[119] yulab.utils_0.0.5         digest_0.6.29            
[121] graph_1.74.0              showtextdb_3.0           
[123] Biostrings_2.64.1         rmarkdown_2.17           
[125] cellranger_1.1.0          fastmatch_1.1-3          
[127] tidytree_0.4.1            GSEABase_1.58.0          
[129] DelayedMatrixStats_1.18.1 nlme_3.1-160             
[131] lifecycle_1.0.3           jsonlite_1.8.2           
[133] Rhdf5lib_1.18.2           BiocNeighbors_1.14.0     
[135] viridisLite_0.4.1         fansi_1.0.3              
[137] pillar_1.8.1              lattice_0.20-45          
[139] KEGGREST_1.36.3           fastmap_1.1.0            
[141] httr_1.4.4                GO.db_3.15.0             
[143] glue_1.6.2                iterators_1.0.14         
[145] png_0.1-7                 bluster_1.6.0            
[147] bit_4.0.4                 HDF5Array_1.24.2         
[149] ggforce_0.4.1             stringi_1.7.8            
[151] sass_0.4.2                blob_1.2.3               
[153] BiocSingular_1.12.0       memoise_2.0.1            
[155] irlba_2.3.5.1             ape_5.6-2                
[157] sysfonts_0.8.8           
</div>
