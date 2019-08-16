//
// Created by Vitor Hadad on 2019-08-15.
//

#include "main.h"
#include <iostream>
#include <stdio.h>
#include <memory>
#include <algorithm>
#include <math.h>
#include <chrono>
#include <ctime>
#include "main.h"


int main(int argc, char **argv) {

    auto time_start = std::chrono::system_clock::now();

    CLI::App app{"grf"};

    std::string forest_type;
    std::string train_filename;
    std::string test_filename;
    std::string oob_output_filename;
    size_t outcome_index = 0;
    size_t treatment_index = 1;
    bool estimate_variance = true;

    // Parsing arguments
    app.add_option("forest_type", forest_type, "Type of forest in {regression, causal}");
    app.add_option("--train", train_filename, "Training data file");
    app.add_option("--test", test_filename, "Test data file");
    app.add_option("--oob_predictions", oob_output_filename, "OOB predictions file");
    app.add_option("--test_predictions", outcome_index, "Test predictions file");
    app.add_option("--outcome_index", outcome_index, "Outcome column number (Zero-indexed)");
    app.add_option("--treatment_index", treatment_index, "Treatment column number (Zero-indexed)");

    CLI11_PARSE(app, argc, argv)

    std::cout << "Fitting forest_type: " << forest_type << "\n";
    std::cout << "Training data set: " << train_filename << "\n";

    train_filename = "/Users/vitorh/Documents/grf/core/test/forest/resources/regression_data.csv";

    // Load data
    std::shared_ptr<Data> train_data(load_data(train_filename));
    std::shared_ptr<Data> test_data;
    if (app.count("--test")) {
        test_data = std::shared_ptr<Data>(load_data(train_filename));
    } else {
        test_data = train_data;
    }


    // Fit and compute oob predict predictions
    if (forest_type == "regression") {
        train_data->set_outcome_index(outcome_index);
        size_t num_covariates = train_data->get_num_cols() - 1;
        ForestOptions options = default_options(num_covariates);
        ForestTrainer trainer = ForestTrainers::regression_trainer();
        Forest forest = trainer.train(train_data.get(), options);
        ForestPredictor predictor = ForestPredictors::regression_predictor(options.get_num_threads());
        std::vector<Prediction> predictions = predictor.predict_oob(forest, train_data.get(), estimate_variance);

        std::cout << std::fixed;
        std::cout << std::setprecision(22);
        for (auto&p : predictions) {
            std::cout << p.get_predictions()[0] << "\t" <<  p.get_variance_estimates()[0] << "\n";
        }

    } else if (forest_type == "causal") {

        size_t num_covariates = train_data->get_num_cols() - 2;

        // Forest options (common for all forests)
        ForestOptions options = default_options(num_covariates);

        // Train outcome forest
        train_data->set_outcome_index(outcome_index);
        train_data->set_treatment_index(treatment_index); // Hack to update disallowed variables
        ForestTrainer outcome_trainer = ForestTrainers::regression_trainer();
        Forest outcome_forest = outcome_trainer.train(train_data.get(), options);
        ForestPredictor outcome_predictor = ForestPredictors::regression_predictor(options.get_num_threads());
        std::vector<Prediction> outcome_predictions = outcome_predictor.predict_oob(outcome_forest, train_data.get(), false);

        // Train treatment forest
        train_data->set_outcome_index(treatment_index);
        train_data->set_treatment_index(outcome_index);  // Hack to update disallowed variables
        ForestTrainer treatment_trainer = ForestTrainers::regression_trainer();
        Forest treatment_forest = treatment_trainer.train(train_data.get(), options);
        ForestPredictor treatment_predictor = ForestPredictors::regression_predictor(options.get_num_threads());
        std::vector<Prediction> treatment_predictions = outcome_predictor.predict_oob(treatment_forest, train_data.get(), false);

        // Replace (W, Y) with (W.hat, Y.hat)
        reset_values(train_data.get(), outcome_index, outcome_predictions);
        reset_values(train_data.get(), treatment_index, treatment_predictions);

        // Train causal forest and compute oob predictions
        train_data->set_outcome_index(outcome_index);
        train_data->set_treatment_index(treatment_index);
        train_data->set_instrument_index(treatment_index);
        ForestTrainer causal_trainer = ForestTrainers::instrumental_trainer(0., false);
        Forest causal_forest = causal_trainer.train(train_data.get(), options);
        ForestPredictor predictor = ForestPredictors::instrumental_predictor(options.get_num_threads());
        std::vector<Prediction> causal_predictions = predictor.predict_oob(causal_forest, train_data.get(), estimate_variance);

        std::cout << std::fixed;
        std::cout << std::setprecision(22);
        for (auto&p : causal_predictions) {
            std::cout << p.get_predictions()[0] << "\t" <<  p.get_variance_estimates()[0] << "\n";
        }

    }

    auto time_end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = time_end-time_start;

    std::cout << "Finished in: " << elapsed_seconds.count() << "s\n";

    return 0;
}


