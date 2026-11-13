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

4. Getis-Ord Gi* computed with 9999 permutations and a fixed seed (42).
A KNN (k=5) spatial weight matrix was used to capture local variance, with binary weights (transform="B").
After FDR correction (alpha=0.2), the significant EDs are:
12 Hotspots (Z > 1.96) and 8 Coldspots (Z < -1.96).

5. Map saved as GiStar_ExtensionDensity_DublinED.png with a diverging colour scheme.

6. Comparison with LISA:
Both methods generally agree on the main clusters. According to the cross-tabulation, out of the 14 High-High areas identified by LISA, 11 are successfully captured as Gi* Hotspots.
Similarly, Gi* captures 7 of the 29 Low-Low areas as Coldspots.
However, there are some differences. Gi* calculates a local sum, which might mask spatial outliers (e.g., High-Low or Low-High), while LISA effectively identifies them.
Additionally, LISA used Queen weights while Gi* used KNN (k=5), which restricts the neighborhood and alters local variance."

"""

from pathlib import Path
import os

# Fix for Windows path issue with non-ASCII usernames (Wrong Chinese characters)
os.environ["JOBLIB_TEMP_FOLDER"] = "C:/temp_joblib"

import numpy as np
import pandas as pd
import geopandas as gpd
import libpysal
import esda
import matplotlib.pyplot as plt
from libpysal.weights import KNN
from esda.getisord import G_Local
from esda import fdr
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
SEED = 42
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

# 1. Use the existing dataset from Q3.1 (which already includes lisa_class)
gdf_gi = eds.copy()

# 2. Apply Log transformation to reduce skewness in the density data
gdf_gi["log_density"] = np.log1p(gdf_gi["extension_density"])
x = gdf_gi["log_density"].values

# 3. Build spatial weights
# Using KNN (k=5) to capture local variance better than Queen contiguity here
w = KNN.from_dataframe(gdf_gi, k=5)
w.transform = "B"

# 4. Compute Getis-Ord Gi*
# Note: 9999 permutations with a fixed seed (42) for reproducibility
print("Computing Gi* (9999 permutations)...")
gi = G_Local(x, w, transform="B", star=1.0, permutations=9999, seed=42)

# 5. Extract Z-scores and calculate two-sided p-values
gdf_gi["gi_z"] = gi.Zs
p_sim = gi.p_sim
p_two_sided = 2 * np.minimum(p_sim, 1 - p_sim)
gdf_gi["gi_p"] = np.clip(p_two_sided, 0, 1)

# 6. Apply False Discovery Rate (FDR) correction
# Alpha set to 0.2 to counteract the extreme strictness of the correction
threshold = fdr(gdf_gi["gi_p"].values, 0.2)
gdf_gi["gi_sig"] = gdf_gi["gi_p"] <= threshold
print(f"FDR Threshold: {threshold:.5f} | Significant EDs: {gdf_gi['gi_sig'].sum()}")


# 7.  Classify EDs into Hotspots, Coldspots, or Not Significant
def classify_gi(row):
    if not row["gi_sig"]:
        return "Not Significant"
    elif row["gi_z"] > 1.96:
        return "Hotspot"
    elif row["gi_z"] < -1.96:
        return "Coldspot"
    return "Not Significant"


gdf_gi["gi_class"] = gdf_gi.apply(classify_gi, axis=1)
print("\nGi* Classification counts:")
print(gdf_gi["gi_class"].value_counts())

# Save the result table for Q3.6 comparison
gdf_gi.drop(columns="geometry").to_csv(OUTPUTS / "GiStar_results.csv", index=False)


#############################################
##### 5. Map the hotspots and coldspots #####
#############################################

# 1. Set up the figure
fig, ax = plt.subplots(figsize=(8, 8))

# 2. Define a diverging colour scheme (Red for Hotspots, Blue for Coldspots, Grey for Not Significant)
colors = {"Hotspot": "#d73027", "Coldspot": "#4575b4", "Not Significant": "#d9d9d9"}

# 3. Plot each class
for cls, color in colors.items():
    sub = gdf_gi[gdf_gi["gi_class"] == cls]
    if not sub.empty:
        sub.plot(ax=ax, color=color, edgecolor="white", linewidth=0.3, label=cls)

# 4. Add title and legend, remove axis for a clean map
ax.set_title("Getis-Ord Gi*: Extension Density Hotspots/Coldspots")
ax.legend(loc="lower left")
ax.set_axis_off()

# 5. Save the map as a PNG file
plt.savefig(
    OUTPUTS / "GiStar_ExtensionDensity_DublinED.png", dpi=200, bbox_inches="tight"
)
print("\nMap saved to outputs/GiStar_ExtensionDensity_DublinED.png")


####################################
##### 6. Compare Gi* with LISA #####
####################################

print("\n--- Q3.6 Cross-tabulation ---")
crosstab_result = pd.crosstab(gdf_gi["lisa_class"], gdf_gi["gi_class"])
print(crosstab_result)

# output the result
crosstab_result.to_csv(OUTPUTS / "LISA_Gi_Comparison.csv")
print("Comparison table saved to outputs/LISA_Gi_Comparison.csv")
