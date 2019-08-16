//
// Created by Vitor Hadad on 2019-08-15.
//

#ifndef GRF_MAIN_H
#define GRF_MAIN_H

#include "commons/utility.h"
#include "forest/ForestOptions.h"
#include "forest/ForestPredictor.h"
#include "forest/ForestPredictors.h"
#include "forest/ForestTrainer.h"
#include "forest/ForestTrainers.h"
#include "CLI/CLI.hpp"

ForestOptions default_options(uint p) {
    bool honesty = true;
    size_t ci_group_size = 2;
    double honesty_fraction = 0.5;
    uint num_trees = 2000;
    double sample_fraction = 0.5;
    uint mtry = std::min(static_cast<uint>(std::ceil(sqrt(p) + 20)), p);
    uint min_node_size = 5;
    double alpha = 0.0;
    double imbalance_penalty = 0.0;
    std::vector<size_t> empty_clusters;
    uint samples_per_cluster = 0;
    uint num_threads = 1;
    uint seed = 42;
    bool prune_empty_leaves = true;

    return ForestOptions(
            num_trees,
            ci_group_size,
            sample_fraction,
            mtry,
            min_node_size,
            honesty,
            honesty_fraction,
            prune_empty_leaves,
            alpha,
            imbalance_penalty,
            num_threads,
            seed,
            empty_clusters,
            samples_per_cluster);
}


void reset_values(Data *data, size_t column_index, std::vector<Prediction> predictions) {
    bool b = false;
    for (auto i = 0; i < predictions.size(); ++i) {
        data->set(column_index, i, predictions[i].get_predictions()[0], b);
    }
}


#endif //GRF_MAIN_H
