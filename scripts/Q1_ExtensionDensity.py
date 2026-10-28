"""
GEOG41160 Assignment 2, Question 1: Extension density by ED
Team: team 1

ANSWERS
Write every paragraph, report, or interpretation the question asks for here,
numbered to match the question's items. Quote the values this script prints,
and update the answers whenever the code changes. For example, for
Question 2:

1. Distribution: skewness of tree density is x.xx, so we used np.log1p because ...
2. Weights: queen contiguity, because ...
3. Moran's I: raw 0.xx (p = 0.xxx); log1p 0.xx (p = 0.xxx),
   9,999 permutations, seed 42.
4. Interpretation: ...
"""

from pathlib import Path

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
## Add the libraries this question needs, for example:
## from libpysal.weights import Queen
## from esda.moran import Moran


#########################################
##### SETUP #####
#########################################

ROOT = Path(__file__).resolve().parents[1]  # the repository folder
DATA = ROOT / "data"
OUTPUTS = ROOT / "outputs"
OUTPUTS.mkdir(exist_ok=True)

SEED = 42  # the same seed in every script


####################################################
##### 1. Create a DuckDB spatial database file #####
####################################################

## What this step does, in one or two lines.

# your code here
# print(...)                                    # print every value the answers quote


###################################################################
##### 2. Print the CRS, number of features, and geometry type #####
###################################################################

# your code here


#################################################################################################
##### 3. Spatially join applications to the EDs in SQL and cout the applications in each ED #####
#################################################################################################

# your code here

########################################################################
##### 4. Calculate the density of applications in each ED and save #####
########################################################################

# your code here

################################################
##### 5. Map the density, and save the map #####
################################################

# your code here

## Save outputs with the file names the brief asks for, for example:
## fig.savefig(OUTPUTS / "LISA_TreeDensity_DublinED.png", dpi=200, bbox_inches="tight")
