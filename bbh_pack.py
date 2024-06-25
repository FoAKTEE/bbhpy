# %load_ext autoreload
# %autoreload 2

import os
import glob
import pathlib
import sys
import argparse
import warnings
import h5py
import numpy as np
import yt

import time
import scipy

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as mpatches
import matplotlib.ticker as mticker
from matplotlib import rcParams
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import multiprocessing
from multiprocessing import Process, Lock, Array
from pathlib import Path

sys.path.append("/Users/hyw/erm/ppscript/vis/python/") # local
sys.path.append("/Users/hyw/Downloads/AthenaK-sink/") # local
sys.path.append("/home/haiyangw/vis/python/") # server: anta

import bin_convert as bc
import athena_read

plt.style.use('bbh.mplstyle') # from C.G. Kim's ncr_paper/tigress mplstyle

# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": ["Computer Modern Sans Serif"]})
