import deepchem as dc
import numpy as np

x = np.random.random((4, 5))
y = np.random.random((4, 1))

dataset = dc.data.NumpyDataset(x, y)

np.array_equal(x, dataset.X)

np.array_equal(y, dataset.Y)

import rdkit

tox21_tasks, tox21_datasets, transformers = dc.molnet.load_tox21()

tox21_tasks

len(tox21_tasks)

tox21_datasets

train_dataset, valid_dataset, test_dataset, = tox21_datasets

train_dataset.X.shape

valid_dataset.X.shape

test_dataset.X.shape

np.shape(train_dataset.Y)

np.shape(valid_dataset.Y)

np.shape(test_dataset.Y)

train_dataset.w.shape

np.count_nonzero(train_dataset.w)

np.count_nonzero(train_dataset.w == 0)

transformers

model = dc.models.MultitaskClassifier(n_tasks=12, n_features=1024, layer_sizes=[1000])

model.fit(train_dataset, nb_epoch=10)

metric = dc.metrics.Metric(dc.metrics.roc_auc_score, np.mean)

import tensorflow as tf

train_scores = model.evaluate(train_dataset, [metric], transformers)
test_scores = model.evaluate(train_dataset, [metric], transformers)

print(train_scores)
print(test_scores)

roc_curve