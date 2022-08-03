# scTEP
devtools::install_github('GRYGY1215/scTEP')
```
data('genesets')

data('goolam')

data = preprocessing(goolam)
data_fa = fa(data, genesets, data_org = 'mmu')
allCluster = clustering(data)
out = trajectoryinference(data, 'mmu', data_fa, allCluster)

scTEP_plot('goolam', data, scDHA_res, out, 'dev_pseudotime', 'umap')
```
