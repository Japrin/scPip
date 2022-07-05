#!/usr/bin/env python3

import argparse
import os
import vireoSNP
import numpy as np
from scipy import sparse
from scipy.io import mmread
import matplotlib.pyplot as plt
import pickle
import pathlib

parser = argparse.ArgumentParser(description='wrapper for running vireoSNP')
parser.add_argument("-n",'--nClone', dest='opt_nClone', type=int,
                    default=3, help='an integer for the accumulator')
parser.add_argument("-i",'--inDir', dest='in_dir', required=True, type=pathlib.Path,
                    help=('directory where the input data (passed_ad.mtx,'
                          'passed_dp.mtx) are stored.'))
parser.add_argument("-o",'--outPrefix', dest='out_prefix', required=True,
                    type=str, help=('output prefix'))

args = parser.parse_args()

print(vireoSNP.__version__)

#print(args)

opt_nClone = args.opt_nClone
in_dir = args.in_dir
out_prefix = args.out_prefix

#opt_nClone = 2
#in_dir = "./t.OUT/mquad/"
#out_prefix = "./t.OUT/mquad/vireo"

os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

AD = mmread("%s/passed_ad.mtx" % in_dir).tocsc()
DP = mmread("%s/passed_dp.mtx" % in_dir).tocsc()

#### model fitting
from vireoSNP import BinomMixtureVB
_model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=opt_nClone)
_model.fit(AD, DP, min_iter=30, n_init=50)
print(_model.ELBO_iters[-1])

#### save result
with open("%s.model.pkl" %  out_prefix, 'wb') as out_file:
    pickle.dump(_model, out_file)

#### visualization
fig = plt.figure(figsize=(11, 4))
plt.subplot(1, 2, 1)
plt.hist(_model.ELBO_inits)
plt.ylabel("Frequency")
plt.xlabel("ELBO in multiple initializations")

plt.subplot(1, 2, 2)
plt.plot(_model.ELBO_iters)
plt.xlabel("Iterations")
plt.ylabel("ELBO in a single initialization")

plt.tight_layout()
#plt.show()
plt.savefig("%s.ELBO.00.pdf" % out_prefix)


# In mitochondrial, allele frequency is highly informative between 0.01 to 0.1,
# so we rescale the colour to give more spectrum for this region.
# You can design/choose your own colors from here:
# https://matplotlib.org/stable/tutorials/colors/colormaps.html

from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

raw_col = cm.get_cmap('pink_r', 200)
new_col = np.vstack((raw_col(np.linspace(0, 0.7, 10)),
                     raw_col(np.linspace(0.7, 1, 90))))
segpink = ListedColormap(new_col, name='segpink')

from vireoSNP.plot import heat_matrix

fig = plt.figure(figsize=(6, 4), dpi=100)
plt.subplot(1, 2, 1)
im = heat_matrix(_model.ID_prob, cmap="Blues", alpha=0.8,
                 display_value=False, row_sort=True)
plt.colorbar(im, fraction=0.046, pad=0.04)
plt.title("Assignment probability")
plt.xlabel("Clone")
plt.ylabel("%d cells" %(_model.n_cell))
plt.xticks(range(_model.n_donor))

plt.subplot(1, 2, 2)
im = heat_matrix(_model.beta_mu, cmap=segpink, alpha=0.8,
                 display_value=False, row_sort=True)
plt.colorbar(im, fraction=0.046, pad=0.04)
plt.title("Mean allelic ratio")
plt.xlabel("Clone")
plt.ylabel("%d SNPs" %(_model.n_var))
plt.xticks(range(_model.n_donor))

plt.tight_layout()
#plt.show()
plt.savefig("%s.prob.allelicRatio.00.pdf" % out_prefix)

#### diagnosis
## repet 5 times, check the ELBO
print("repeat model fitting 5 times, see whether the same (best) ELBO is found")
n_init = 50
for i in range(5):
    _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=opt_nClone)
    _model.fit(AD, DP, min_iter=30, n_init=n_init)
    print("rerun %d:" %i, _model.ELBO_iters[-1])

## check whether the number of clones choosed is good
print("check whether the number of clones choosed is good")
n_init = 50
n_clone_list = np.arange(2, 6)

_ELBO_mat = []
for k in n_clone_list:
    _model = BinomMixtureVB(n_var=AD.shape[0], n_cell=AD.shape[1], n_donor=k)
    _model.fit(AD, DP, min_iter=30, n_init=n_init)
    _ELBO_mat.append(_model.ELBO_inits)

fig = plt.figure(figsize=(6, 4))
plt.plot(np.arange(1, len(n_clone_list)+1), np.max(_ELBO_mat, axis=1))
plt.boxplot(_ELBO_mat)
plt.xticks(np.arange(1, len(n_clone_list)+1), n_clone_list)
plt.ylabel("ELBO")
plt.xlabel("n_clones")
plt.savefig("%s.ELBO.nClone.00.pdf" % out_prefix)
#plt.show()

##############
#mtSNP_ids = ['mt_variant%d' %x for x in range(AD.shape[0])]
#cell_label = np.array(['clone1'] * 27 + ['clone2'] * 27 + ['clone3'] * 27)
#id_uniq = ['clone1', 'clone2', 'clone3']
#vireoSNP.plot.anno_heat(AD/DP, col_anno=cell_label, col_order_ids=id_uniq,
#                        cmap=segpink, yticklabels=mtSNP_ids)


