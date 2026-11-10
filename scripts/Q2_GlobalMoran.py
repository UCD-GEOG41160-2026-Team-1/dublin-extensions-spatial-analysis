"""
GEOG41160 Assignment 2, Question 2: Global spatial autocorrelation
Team: team 1

ANSWERS

1. Distribution:
The original extension density was positively skewed (1.391).
A log1p transformation produced skewness of -1.140,
while a square-root transformation reduced skewness
to 0.354. The square-root transformation was selected because
it produced more symetrical distribution.

2. Spatial weights:
Queen contiguity weights were used because EDs sharing
a boundary or vertex are considered neighbours.
The analysis included 162 EDs, with an average of
5.70 neighbours per ED. There were no islands and
one connected component. Weights were row-standardised.

3. Global Moran's I:
Using 999 permutations and random seed 42:

Raw extension density:
Moran's I = 0.332, p = 0.001.

Square-root extension density:
Moran's I = 0.375, p = 0.001.

Expected Moran's I = -0.006.

4. Interpretation:
Both variables exhibited statistically significant
positive spatial autocorrelation (p = 0.001).
This indicates that EDs with similar extension
densities tend to occur near one another.

The square-root transformation increased Moran's I
from 0.332 to 0.375, indicating a slightly stronger
measured spatial association. The transformation
reduced the influence of extreme density values but
did not change the conclusion that the spatial
pattern is clustered.
"""

from pathlib import Path

import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
## Add the libraries this question needs, for example:
from libpysal.weights import Queen, lag_spatial
from esda.moran import Moran


#########################################
##### SETUP #####
#########################################

ROOT = Path(__file__).resolve().parents[1]  # the repository folder
DATA = ROOT / "data"
OUTPUTS = ROOT / "outputs"
OUTPUTS.mkdir(exist_ok=True)

SEED = 42

EXTENSION_DENSITY = OUTPUTS / "DublinED_ExtensionDensity.gpkg"
if not EXTENSION_DENSITY. exists():
    raise FileNotFoundError(
        f"{EXTENSION_DENSITY} not found: run Q1_ExtensionDensity.py first."
    )
eds = gpd.read_file(EXTENSION_DENSITY)

print(len(eds), "EDs,", eds.crs)

##########################################################################
##################1.Distribution of extension density#######################
##########################################################################

density=eds["extension_density"].astype(float)
print("Valid EDs:", density.notna().sum())
print("Missing values:",density.isna().sum())
print("Minimum density:", density.min())
print("Maximum density:", density.max())
print("Mean density:", density.mean())
print("Skewness:", density.skew())

fig, ax=plt.subplots(figsize=(8,5))

ax.hist(density.dropna(), bins=20, edgecolor="black")
ax.set_title("Distribution of Extension Density in Dublin City")
ax.set_xlabel("Extension Applications per km2")
ax.set_ylabel("Number of Electoral Divisions")

fig.tight_layout()

fig.savefig(
    OUTPUTS / "ExtensionDensity_Histogram.png",
    dpi=200,
    bbox_inches="tight"
)
plt.close(fig)
print("Histogram saved successfully.")

density_log=np.log1p(density)

print("Original skewness:", density.skew())
print("Log1p skewness:", density_log.skew())

density_sqrt=np.sqrt(density)
print("Square-root skewness:",density_sqrt.skew())

# Compare original and square-root transformed extension density

fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# Original distribution
axes[0].hist(density, bins=30, color="#4a7c3a")
axes[0].set_title("Original Extension Density")
axes[0].set_xlabel("Applications per km²")
axes[0].set_ylabel("Number of EDs")

# Square-root transformed distribution
axes[1].hist(density_sqrt, bins=30, color="#4a7c3a")
axes[1].set_title("Square-root Extension Density")
axes[1].set_xlabel("Square root of applications per km²")
axes[1].set_ylabel("Number of EDs")

plt.tight_layout()

fig.savefig(
    OUTPUTS / "ExtensionDensity_Comparison.png",
    dpi=200,
    bbox_inches="tight"
)

plt.close(fig)

print("Comparison histogram saved successfully.")

# Part 1, Step 3: Visualise transformed extension density

eds["sqrt_density"] = density_sqrt

fig, ax = plt.subplots(figsize=(7, 6))

eds.plot(
    column="sqrt_density",
    cmap="YlGn",
    legend=True,
    edgecolor="grey",
    linewidth=0.3,
    ax=ax,
    legend_kwds={"label": "Square root of applications per km²"}
)

ax.set_title("Extension Density by ED, Dublin City (Square-root Scale)")
ax.set_axis_off()

plt.tight_layout()

fig.savefig(
    OUTPUTS / "ExtensionDensity_Sqrt_Map.png",
    dpi=180,
    bbox_inches="tight"
)

plt.close(fig)

print("Transformed density map saved successfully.")


#############################################################################
##### 2. Build spatial weights, and state which definition of neighbour #####
#############################################################################
w = Queen.from_dataframe(eds, use_index=False)

# Check the spatial weights
print("Number of EDs in weights:", w.n)
print("Average neighbours per ED:", w.mean_neighbors)
print("EDs without neighbours:", w.islands)
print("Connected components:", w.n_components)

# Row-standardise the spatial weights
w.transform = "r"
# Part 1, Step 5: Calculate spatial lag

# Use square-root transformed extension density
x = density_sqrt.to_numpy()

# Calculate the average density of neighbouring EDs
lag = lag_spatial(w, x)

print("Citywide average:", round(x.mean(), 3))
print("Number of spatial lag values:", len(lag))
print("First five ED values:", np.round(x[:5], 3))
print("First five neighbour averages:", np.round(lag[:5], 3))

print("Weights transformation:", w.transform)



############################################################################
#3: Calculate Global Moran's I##############################################
############################################################################
# Set the random seed for reproducibility
np.random.seed(SEED)

# Calculate Moran's I using square-root density
moran_sqrt = Moran(x, w, permutations=999)

# Print results
print("Square-root Moran's I:", round(moran_sqrt.I, 3))
print("Expected Moran's I:", round(moran_sqrt.EI, 3))
# Part 3: Moran Scatterplot

# Standardise the square-root density values
z = (x - x.mean()) / x.std()

# Calculate spatial lag of standardised values
z_lag = lag_spatial(w, z)

# Create Moran scatterplot
fig, ax = plt.subplots(figsize=(5.5, 5))

ax.scatter(z, z_lag, s=12)

# Calculate and plot regression slope
slope = np.polyfit(z, z_lag, 1)[0]

ax.plot(
    [z.min(), z.max()],
    [slope * z.min(), slope * z.max()],
    color="red"
)

# Add reference lines
ax.axhline(0, color="grey")
ax.axvline(0, color="grey")

# Label the four quadrants
corners = {
    "HH": (0.97, 0.97),
    "LH": (0.03, 0.97),
    "LL": (0.03, 0.03),
    "HL": (0.97, 0.03)
}

for label, (cx, cy) in corners.items():
    ax.text(
        cx, cy, label,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        color="grey",
        ha="right" if cx > 0.5 else "left",
        va="top" if cy > 0.5 else "bottom"
    )

ax.set_xlabel("ED's own value (standardised)")
ax.set_ylabel("Average of its neighbours")
ax.set_title(f"Moran scatterplot: slope = {slope:.3f}")

plt.tight_layout()

fig.savefig(
    OUTPUTS / "Moran_Scatterplot_ExtensionDensity.png",
    dpi=200,
    bbox_inches="tight"
)

plt.close(fig)

print("Moran scatterplot saved successfully.")
print("Moran scatterplot slope:", round(slope, 3))
# Part 4: Significance of Global Moran's I

# Print permutation test results
print("Permutation p-value:", round(moran_sqrt.p_sim, 3))

print(
    "Random maps at least as clustered:",
    (moran_sqrt.sim >= moran_sqrt.I).sum(),
    "of",
    len(moran_sqrt.sim)
)

# Visualise the permutation distribution
fig, ax = plt.subplots(figsize=(6, 3.5))

ax.hist(moran_sqrt.sim, bins=30, color="grey")

ax.axvline(
    moran_sqrt.I,
    color="red",
    label="Observed Moran's I"
)

ax.set_xlabel("Moran's I")
ax.set_ylabel("Frequency")
ax.set_title("999 Random Permutations vs Observed Moran's I")

ax.legend()

plt.tight_layout()

fig.savefig(
    OUTPUTS / "Moran_Permutation_Test.png",
    dpi=200,
    bbox_inches="tight"
)

plt.close(fig)

print("Permutation test histogram saved successfully.")
# Compare raw and transformed extension density

# Use the same random seed and spatial weights
np.random.seed(SEED)

# Calculate Moran's I for original extension density
moran_raw = Moran(density.to_numpy(), w, permutations=999)

# Print results for comparison
print("Raw density Moran's I:", round(moran_raw.I, 3))
print("Raw density p-value:", round(moran_raw.p_sim, 3))

print("Square-root Moran's I:", round(moran_sqrt.I, 3))
print("Square-root p-value:", round(moran_sqrt.p_sim, 3))