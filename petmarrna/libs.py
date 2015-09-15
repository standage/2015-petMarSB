#!/usr/bin/env ipython
from __future__ import print_function

# Standard library packages

from itertools import izip
import os
import sys
import json
from collections import defaultdict as ddict
from pprint import pprint
import warnings

# Dependencies

from bioservices import UniProt

from IPython.display import SVG
from IPython.display import FileLink
from IPython.display import HTML

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles

import mygene

import numpy as np
from numpy.random import rand

import pandas as pd
from pandas import Categorical

import screed

from scipy import stats
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch

import seaborn as sns

from sklearn.decomposition import PCA, KernelPCA
from sklearn.cluster import KMeans

# Local modules

import blasttools
import buscotools
import hmmertools
from figuremanager import FigManager


