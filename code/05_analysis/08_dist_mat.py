import numpy as np
from sklearn.metrics import pairwise_distances_chunked
import pandas as pd
df = pd.read_csv("data/human/04_sl_input/cell_line_compare_expr.csv")

X = df.iloc[:,1:10721].transpose()
gen = pairwise_distances_chunked(X)
gen2 = next(gen) # this was enough
#gen3 = next(gen) 

np.savetxt("cell_dist_mat.txt", gen2)
