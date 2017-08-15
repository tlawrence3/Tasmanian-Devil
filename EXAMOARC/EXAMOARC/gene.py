from __future__ import print_function
import pandas as pd
import numpy as np

def gene_classify(file_name, upper, lower, file_out):
    geneDF = pd.read_csv(file_name, sep = "\t", header=None, names = ["Gene", "Count"])
    geneDF["class"] = np.where(geneDF["Count"] >= geneDF.Count.quantile(upper),
                               1, np.where(geneDF["Count"] <= geneDF.Count.quantile(lower), -1, 0))
    geneDF.to_csv(file_out, sep = "\t", columns=["Gene","class"], header = False, index = False)

if __name__ == "__main__":
    gene_classify("test.txt", .75, .25, "test.counts.txt")

