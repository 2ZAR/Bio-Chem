#Practice#
#ECFP#

smiles = ['C1CCCC1', 'O1CCOCC1']
mols = [Chem.MolFromSmiles(smile) for smile in smiles]
feat = dc.feat.CircularFingerprint(size=1024)
arr = feat.featurize(mols)

feat = dc.feat.RDKitDescriptors()
arr = feat.featurize(mols)

tasks, datasets, transformers = dc.Molnet.load_delaney(featurizer='GraphConv')
train_dataset, valid_dataset, test_dataset = datasets

model = GraphConvModel(n_tasks=1, mode='regression', dropout=0.2)
model.fit(train_dataset, np_epoch=100)

metric = dc.metrics.Metric(dc.metrics.pearson_r2_score)
print(model.evaluate(train_dataset, [metric], transformers))
print(model.evaluate(test_dataset, [metric], transformers))

smiles = [
    'COC(C)(C)CCCC(C)CC=CC(C)=CC(=0)OC(C)C',
    'CCOC(=0)CC'
    'CSc1nc(NC(C)C)nc(NC(C)C)n1'
    'CC(C#C)N(C)C(=0)Nc1ccc(Cl)cc1'
    'Cc1cc2ccccc2cc1C']

from rdkit import Chem
mols = [Chem.MolFromSmiles(s) for s in smiles]
featurizer = dc.feat.ConvMolFeaturizer()
x = featurizer.featurize(mols)

predicted_solubility = model.predict_on_batch(x)

