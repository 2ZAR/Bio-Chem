#Practice#
#PDBBind dataset featurize - RandomForest#

import deepchem as dc

grid_featurizer = dc.feat.RdkitGridFeaturizer(voxel_width=2.0,
                                              feature_types =['hbond', 'salt_bridge', 'pi_stack', 'cation_pi', 'ecfp', 'splif']
                                              sanitize=True, flatten=True)

tasks, datasets, transformers = dc.molnet.load_pdbbind(featurizer="grid", split="random", subset="core")
train_dataset, valid_dataset, test_dataset = datasets

import sklearn as sk
from sklearn.ensemble import RandomForestRegressor

sklearn_model = RandomForestRegressor(n_estimators=100)
model = dc.models.SklearnModel(sklearn_model)
model.fit(train_dataset)

n_features = train_dataset.X.shape[1]
model = dc.models.MultitaskRegressor(n_tasks=len(pdbbind_tasks),
                                     n_features=n_features,
                                     layers_size=[2000, 1000],
                                     dropouts=0.5
                                     learning_rate=0.0003)

model.fit(train_dataset, nb_epoch=250)

metric = dc.metrics.Metric(dc.metrics.pearson_r2_score)

print("Evaluating model")

train_scores = model.evaluate(train_dataset, [metric], transformers)
test_scores = model.evaluate(test_dataset, [metric], transformers)

print("Train scores")
print(train_scores)
print("Test scores")
print(test_scores)

