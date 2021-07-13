#!/usr/bin/env Rscript

# Load optparse we need to check inputs
suppressPackageStartupMessages(require(optparse))
# Load common functions
suppressPackageStartupMessages(require(workflowscriptscommon))

option_list = list(
    make_option(
        c("-i", "--input-classifier-object"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the input classifier object in .rds format.'
  ),
    make_option(
        c("-q", "--query-expression-data"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the SCE object containing expression data to be predicted.'
  ),
    make_option(
        c("-n", "--n-rand-prof"), 
        action = "store",
        default = 2,
        type = 'numeric',
        help = 'The number of random profiles generated for evaluation process'
  ),
    make_option(
        c("-o", "--prediction-output"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path to the predictions obtained from the classifier.'
  ),
    make_option(
        c("-r", "--return-raw-output"), 
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = "Should the output be returned in raw format (i.e. not transformed into table)? Default: FALSE."
  )
)

# parse args 
opt <- wsc_parse_args(option_list, mandatory = c('input_classifier_object', 'query_expression_data', 'prediction_output'))

# load remaining packages 
suppressPackageStartupMessages(require(singleCellNet))
suppressPackageStartupMessages(require(SingleCellExperiment))

# get query matrix
query = readRDS(opt$query_expression_data)
if("counts" %in% names(assays(query))){
    pred_data = counts(query)
} else if("normcounts" %in% names(assays(query))){
    pred_data = normcounts(query)
} else {
    stop("Stopping: Neither 'counts' nor 'normcounts' slot found in provided query object.")
}
classifier = readRDS(opt$input_classifier_object)

# check overlap across clasification genes and query genes; impute mean values if necessary
cl_genes = classifier[["cnProc"]][["cgenes"]] 
query_genes = row.names(pred_data)

if(! all(cl_genes %in% query_genes)) {
    print("...Imputing missing genes")
    missing_genes = cl_genes[which(!cl_genes %in% query_genes)]
    imp_vals = as.numeric(apply(pred_data, 2, mean))
    imp_vals = t(replicate(length(missing_genes), imp_vals))
    row.names(imp_vals) = missing_genes
    pred_data = rbind(pred_data, imp_vals)
}

res = scn_predict(cnProc = classifier[['cnProc']], expDat = pred_data, nrand = opt$n_rand_prof)

if(opt$return_raw_output){
    saveRDS(res, opt$prediction_output)
} else{
    # process output to get a standard table
    res = res[, -grep("rand", colnames(res))]
    cell_id = colnames(res)
    cell_types = row.names(res)
    # get top predictions for each cell 
    .get_top_cells = function(col){
        max_val = max(col)
        max_idx = which.max(col)
        return(c(max_val, max_idx))
    }
    top_cells = apply(res, 2, .get_top_cells)
    scores = as.numeric(top_cells[1,])
    idx = as.integer(top_cells[2,])
    pred_labs = cell_types[idx]
    # build a table 
    tbl = data.frame(cbind(cell_id = cell_id, predicted_label = pred_labs, score = scores))
    dataset = classifier[["dataset"]]
    # add metadata lines
    fc = file(opt$prediction_output)
    writeLines(c("# tool singleCellNet", paste("# dataset", dataset)), fc)
    close(fc)
    write.table(tbl, file = opt$prediction_output, sep="\t", row.names=FALSE, append=TRUE)
}
