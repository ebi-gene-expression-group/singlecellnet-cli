#!/usr/bin/env Rscript

# Load optparse we need to check inputs
suppressPackageStartupMessages(require(optparse))
# Load common functions
suppressPackageStartupMessages(require(workflowscriptscommon))

option_list = list(
    make_option(
        c("-i", "--input-object"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the input object with training data. SCE class in .rds format.'
  ),
    make_option(
        c("-c", "--cell-type-col"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Name of the cell type annotation column in object metadata.'
  ),
    make_option(
        c("-b", "--cell-barcode-col"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Name of the barcode column in object metadata.'
  ),
    make_option(
        c("-n", "--n-top-genes"), 
        action = "store",
        default = 10,
        type = 'numeric',
        help = 'The number of classification genes per category.'
  ),
    make_option(
        c("-p", "--n-top-gene-pairs"), 
        action = "store",
        default = 25,
        type = 'numeric',
        help = 'The number of top gene pairs per category.'
  ),
    make_option(
        c("-r", "--n-rand"), 
        action = "store",
        default = 70,
        type = 'numeric',
        help = 'Number of random profiles to generate for training.'
  ),
    make_option(
        c("-t", "--n-trees"), 
        action = "store",
        default = 1000,
        type = 'numeric',
        help = 'Number of trees for random forest classifier'
  ),
  make_option(
        c("-s", "--stratify"), 
        action = "store_true",
        default = FALSE,
        type = 'logical',
        help = "Should the 'stratify' parameter be set?"
  ),
  make_option(
        c("-e", "--weighted-down-threshold"), 
        action = "store",
        default = 0.25,
        type = 'numeric',
        help = 'The threshold at which anything lower than that is 0 for weighted_down function'
  ),
  make_option(
        c("-f", "--transprop-factor"), 
        action = "store",
        default = 10000,
        type = 'numeric',
        help = "Scaling factor for transprop."
  ),
  make_option(
        c("-w", "--weighted-down-total"), 
        action = "store",
        default = 1500,
        type = 'numeric',
        help = 'Numeric post transformation sum of read counts for weighted_down function'
  ), 
  make_option(
        c("-d", "--dataset-id"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Name of the dataset used for training the classifier.'
  ),
  make_option(
        c("-o", "--output-path"), 
        action = "store",
        default = NA,
        type = 'character',
        help = 'Output path for the trained classifier.'
  )
)

# parse args 
opt <- wsc_parse_args(option_list, mandatory = c('input_object', 'output_path', 'cell_type_col', 'cell_barcode_col'))
# load remaining packages 
suppressPackageStartupMessages(require(singleCellNet))
suppressPackageStartupMessages(require(SingleCellExperiment))

sce_obj = readRDS(opt$input_object)
cell_lab_field = opt$cell_type_col
barcode_field = opt$cell_barcode_col

# subset expression data to remove unlabelled cells 
cells_to_keep = which(! colData(sce_obj)[[cell_lab_field]] %in% c("", NA, "unknown",
                                                                  "rand", "unassigned",
                                                                  "ambiguous")) 
sce_obj = sce_obj[, cells_to_keep]

# We need counts for dropout detection. If normcounts is present, just reassign those
assay_names <- names(assays(sce_obj))
if (! 'counts' %in% assay_names){
    if ( 'normcounts' %in% assay_names){
        names(assays(sce_obj))[names(assays(sce_obj)) == 'normcounts'] <- 'counts'
        assay_names <- names(assays(sce_obj))
    }else{
        stop("Neither 'counts' nor 'normcounts' are populated in input object.")
    }
}

# separate expression matrix and metadata
exp_mat = counts(sce_obj)
exp_metadata = colData(sce_obj)

# run model training 
classifier = scn_train(stTrain = exp_metadata, expTrain = exp_mat, 
                       nTopGenes = opt$n_top_genes, nRand = opt$n_rand, nTrees = opt$n_trees,
                       nTopGenePairs = opt$n_top_gene_pairs, dLevel = opt$cell_type_col,
                       colName_samp = opt$cell_barcode_col, weightedDown_dThresh = opt$weighted_down_threshold,
                       weightedDown_total = opt$weighted_down_total, stratify = opt$stratify)

# add dataset ID, if specified 
if(!is.na(opt$dataset_id)){
    classifier["dataset"] = opt$dataset_id
} else{
    classifier["dataset"] = NA 
}

# save output 
saveRDS(classifier, opt$output_path)
