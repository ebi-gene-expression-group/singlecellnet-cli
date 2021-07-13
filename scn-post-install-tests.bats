#!/usr/bin/env bats

@test "Train a classifier" {
    run rm -rf $classifier && scn-train-model.R\
                                --input-object $training_data\
                                --cell-type-col $cell_type_col\
                                --cell-barcode-col $cell_barcode_col\
                                --output-path $classifier

    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$classifier" ]
}

@test "Predict undersampled data" {
    run rm -rf $prediction_output && scn-predict.R\
                                --input-classifier-object $classifier\
                                --query-expression-data $undersampled_pred_sce\
                                --prediction-output $prediction_output
    
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$classifier" ]


}