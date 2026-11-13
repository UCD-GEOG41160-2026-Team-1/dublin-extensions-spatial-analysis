"""
GEOG41160 Assignment 2, Question 3: Local clusters and hotspots
Team: team 1

ANSWERS
1. The Local Moran's I values were calculated for each ED to identify their local
spatial autocorrelation. Each ED is then classified into one of the 4 quadrants
based on the neighbors of the ED. The quadrants are: High-High, Low-High, Low-Low, and High-Low.
The results show the distribution of EDs across the four quadrants,
with 47 significant results.14 are classified as High-High, 5 as Low-High,
27 as Low-Low, and 1 as High-Low. The remaining 115 EDs are not significant.

2. The false discovery rate correction was applied to the p-values of the Local
Moran's I to minise the amout of false postives/significant results.
No more than 5% was accepted for this work.
30 EDs remain significant after the correction.

3. The LISA classes were mapped, with a legend indicating the number of EDs in
each class. This map allow us to visualise the clusters of extension density
accross the Dublin EDs. As we can see, there is a cluster in the south, with many
EDs classified as High-High, which means a high extension density ED is surrounded
by other EDs with high extension density.The west and the North are concentrating
the Low-Low clusters, which means they don't have high extension density in these EDs.
The map was saved as LISA_ExtensionDensity_DublinED.png in the output folder.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import libpysal
import esda
from esda.moran import Moran_Local

#########################################
##### SETUP #####
#########################################

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUTS = ROOT / "outputs"
OUTPUTS.mkdir(exist_ok=True)


DENSITY_FILE = OUTPUTS / "DublinED_ExtensionDensity.gpkg"
SEED = 42

# This question runs from the question1 outputs, so this checks that the file exists before proceeding.
if not DENSITY_FILE.exists():
    raise FileNotFoundError(
        f"{DENSITY_FILE} not found: run Q1_ExtensionDensity.py first."
    )
eds = gpd.read_file(DENSITY_FILE)

# Inspect dataset to know the datset before analysis
print(len(eds), "EDs,", eds.crs)
print(eds.head())
print(eds.columns)


#############################################
##### 1. Compute Local Moran's I (LISA) #####
#############################################

# Creation of x to hold the extension density values as an array
x = eds["extension_density"].values

# Application of the Queen contiguity weight based on the ED boundaries.
# EDs are neiughbours if they share a common edge or corner.
# The results is stored in the w variable.
w = libpysal.weights.Queen.from_dataframe(eds, use_index=True)

# The weights are standardized, each ED's weights sum to 1. This is applied to the w dataset.
w.transform = "r"

# The Local Moran's I values are then calculated for each ED using the extension
# density values and the spatial weights.
# A random seed is used to ensure the reproducibility of the results.
lisa = esda.Moran_Local(x, w, permutations=9999, seed=SEED)
# It does 9999 permutations to calculate the p-values for each ED's Local Moran's I
#  value, which indicates the significance of the clustering.


##Classify the local Moran's I values into four quadrants based on the mean of the
#  extension density values and significant clusters are identified based on the
# p-values of the Local Moran's I.
quadrants = {1: "High-High", 2: "Low-High", 3: "Low-Low", 4: "High-Low"}
# turns each EDs qudrant into a name
classes = pd.Series(lisa.q).map(quadrants)

# keeps the name when significant, otherwise "Not significant"
eds["lisa_class"] = np.where(lisa.p_sim <= 0.05, classes, "Not significant")
print(eds["lisa_class"].value_counts())


######################################################
##### 2. Apply a false discovery rate correction #####
######################################################

# Apply the false discovery rate correction to prevent false positives

# Allows only 5% of the false discovery to pass through the test
threshold = esda.fdr(lisa.p_sim, 0.05)
print("FDR threshold:", round(threshold, 5), "| smallest p:", round(lisa.p_sim.min()))
print("Significant after FDR correction:", (lisa.p_sim <= threshold).sum())


#####################################################################
##### 3. Map the LISA classes with a legend and save the map ########
#####################################################################
from matplotlib.patches import Patch


# Add a column to the eds with the EDs local values and quadrants classes names
eds["quadrant"] = classes.values

# Define the colours for each class of the map
colours = {
    "High-High": "#d7191c",
    "Low-High": "#abd9e9",
    "Low-Low": "#2c7bb6",
    "High-Low": "#fdae61",
    "Not significant": "#eeeeee",
}

# draw the map
fig, ax = plt.subplots(figsize=(10, 7.5))

eds.plot(color=eds["quadrant"].map(colours), edgecolor="grey", linewidth=0.3, ax=ax)

ax.legend(
    handles=[
        Patch(color=c, label=f"{k} ({(eds['quadrant'] == k).sum()})")
        for k, c in colours.items()
    ],
    loc="lower right",
    fontsize=7,
)

ax.set_title("LISA classes")

# Save the map to the outputs folder
plt.savefig(
    OUTPUTS / " LISA_ExtensionDensity_DublinED.png", dpi=180, bbox_inches="tight"
)
# Show the map
plt.show()


#####################################################################
##### 4. Compute Getis-Ord Gi*, classify hotspots and coldspots #####
#####################################################################

# your code here

#############################################
##### 5. Map the hotspots and coldspots #####
#############################################

# your code here

####################################
##### 6. Compare Gi* with LISA #####
####################################

# your code here

## Save outputs with the file names the brief asks for, for example:
## fig.savefig(OUTPUTS / "LISA_TreeDensity_DublinED.png", dpi=200, bbox_inches="tight")
