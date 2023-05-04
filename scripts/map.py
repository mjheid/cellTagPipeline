import scanpy as sc

ad = sc.read("test.csv", first_column_names=1)

sc.pp.neighbors(ad)
sc.tl.umap(ad)
sc.pl.umap(ad, save="test_umap.csv")