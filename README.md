# singlecellnet-cli
Command Line Interface scripts for the [singleCellNet](https://github.com/pcahan1/singleCellNet) tool. 

## Installation


## Commands 
### scn-train-model
Train a random forest classifier for provided single-cell dataset. Data is accepted as a SingleCellExperiment object. 
```
scn-train-model.R\
    --input-object <Input SCE object in .rds format>\
    --cell-type-col <Name of the cell type annotation column in object metadata>\
    --cell-barcode-col <Name of the barcode column in object metadata>\
    --n-top-genes <The number of classification genes per category>\
    --n-top-gene-pairs <The number of top gene pairs per category>\
    --n-rand <Number of random profiles to generate for training>\
    --n-trees <Number of trees for random forest classifier>\
    --stratify <Should the 'stratify' parameter be set?>\
    --weighted-down-threshold <The threshold at which anything lower than that is 0 for weighted_down function>\
    --transprop-factor <Scaling factor for transprop>\
    --weighted-down-total <Numeric post transformation sum of read counts for weighted_down function>\
    --dataset-id <Name of the dataset used for training the classifier>\
    --output-path <Output path for the trained classifier>
```

### scn-predict
Classify query data using a trained classifier. In case when not all feature genes are present in the query, mean imputation is used. 

```
scn-predict.R\
    --input-classifier-object <Path to the input classifier object in .rds format>\
    --query-expression-data <Path to the SCE object containing expression data to be predicted>\
    --n-rand-prof <The number of random profiles generated for evaluation process>\
    --prediction-output <Output path to the predictions obtained from the classifier>\
    --return-raw-output <Should the output be returned in raw format (i.e. not transformed into table)? Default: FALSE>
```

