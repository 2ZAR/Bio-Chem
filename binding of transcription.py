#Practice#
#CNN : predicting the binding of transcription factors#

model = dc.models.TensorGraph(batch_size=1000)
features = layers.Feature(shape=(None, 101, 4))
labels = layers.Label(shape=(None, 1))
weights = layers.Weights(shape=(None, 1))

prev = features
for i in range(3):
    prev = layers.Conv1D(filters=15, kernel_size=10, activation=tf.nn.relu,
                         padding='same', in_layers=prev)
prev = layers.Dropout(dropout_prob=0.5, in_layers=prev)

logits = layers.Dense(out_channels=1, in_layers=layers.Flatten(prev))
output = layers.Sigmoid(logits)
model.add_output(output)

loss = layers.SigmoidCrossEntropy(in_layers=[labels, logits])
weighted_loss = layers.WeightedError(in_layers=[loss, weights])
model.set_loss(weighted_loss)

train = dc.data.DiskDataset('train_dataset')
valid = dc.data.DiskDataset('valid_dataset')
metric = dc.metrics.Metric(dc.metrics.roc_auc_score)

for i in range(20):
    model.fit(train, nb_epoch=10)
    print(model.evaluate(train, [metric]))
    print(model.evaluate(valid, [metric]))

#Add chromatin accessibility information#

span_accessibility = {}
for line in open('accessibility.txt') :
    fields = line.split()
    span_accessibility[fields[0]] = format(fields[1])

accessibility = layers.Feature(shape=(None, 1))

logits = layers.Dense(out_channels=1, in_layers=layers.Flatten(prev))

prev = layers.Concat([layers.Flatten(prev), accessibility])

logits = layers.Dense(out_channels=1, in_layers=prev)

def generate_batches(dataset, epochs) :
    for epoch in range(epochs) :
        for X, y, w, ids in dataset.iterbatches(batch_size=1000, pad_batches=True) :
            yield {features : X, accessibility : np.array([span_accessibility[id] for id in ids]),
                   labels : y, weights : w}

for i in range(20) :
    model.fit_generator(generate_batches(train, 10))
    print(model.evaluate_generator(generate_batches(train, 1), [metric], labels=[labels], weights=[weights])
    print(model.evaluate_generator(generate_batches(valid, 1), [metric], labels=[labels], weights=[weights])

#RNA interference#

model = dc.models.tensorGraph()
features = layers.Feature(shape=(None, 21, 4))
labels = layers.Label(shape=(None, 1))
prev = features

for i in range(2) :
    prev = layers.Conv1D(filters=10, kernel_size=10, activation=tf.nn.relu, padding='same', in_layers=prev)
    prev = layers.Dropout(dropout_prob=0.3, in_layers=prev)

output = layers.Dense(out_channels=1, activation_fn=tf.sigmoid, in_layers=layers.Flatten(prev))
model.add_output(output)
loss = layers.ReduceMean(layers.L2Loss(in_layers=[labels, output]))
model.set_loss(loss)

train = dc.data.DiskDataset('train_siRNA')
valid = dc.data.DiskDataset('valid_siRNA')
metric = dc.metrics.Metric(dc.metrics.pearsonr, mode='regression')
for i in range(20) :
    model.fit(train, nb_epoch=10)
    print(model.evaluate(train, [metric]))
    print(model.evaluate(valid, [metric]))

    
