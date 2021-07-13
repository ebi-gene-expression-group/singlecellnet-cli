#!/usr/bin/env bash

export PATH=$(pwd):$PATH
test_dir="$(pwd)/post_install_tests"
[ -d $test_dir ] && rm -r $test_dir
mkdir -p $test_dir && pushd $test_dir

# import test data  
wget http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/github_test_data/cli_tools/singlecellnet-cli.tar.gz
tar xvf singlecellnet-cli.tar.gz

export training_data="singlecellnet-cli/SCE.rds"
export classifier="trained_model.rds"
export undersampled_pred_sce="singlecellnet-cli/undersampled_train_sce.rds"
export cell_type_col="Factor.Value..inferred.cell.type...ontology.labels."
export cell_barcode_col="Barcode"
export prediction_output="pred_labels.tsv"

# run tests
scn-post-install-tests.bats 
