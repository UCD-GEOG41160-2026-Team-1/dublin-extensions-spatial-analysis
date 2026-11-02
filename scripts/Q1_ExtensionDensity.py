"""
GEOG41160 Assignment 2, Question 1: Extension density by ED
Team: 1

ANSWERS
1. The national planning table has 510,390 rows. Filtering to Dublin City
   Council leaves 24,879 rows. Keeping applications received in the complete
   calendar years 2019-2025 leaves 23,176 rows. The case-insensitive search
   for "extension" in DevelopmentDescription leaves 9,914 applications.
   The national ED table has 3,420 rows and the Dublin City filter leaves 162.

2. The ED layer has 162 MULTIPOLYGON features in EPSG:2157. The filtered
   application layer has 9,914 POINT features in EPSG:3857 and is reprojected
   to EPSG:2157. Irish Transverse Mercator is a projected CRS in metres, so it
   is suitable for measuring ED areas in Ireland. Web Mercator distorts area
   and distance, which would make the densities inaccurate.

3. The spatial join counts 9,900 applications inside Dublin City EDs. The
   other 14 applications fall outside every ED. Their addresses include places
   outside the city, such as Swords and Skerries. Other points are close to,
   but outside, an ED boundary, so the recorded point locations do not all
   align with the city ED boundaries.

4. Density is the number of applications divided by ED area in square
   kilometres. The output contains all 162 EDs. The left join would keep an
   ED with a count of zero if it had no applications. The result is saved as
   DublinED_ExtensionDensity.gpkg.

5. The density map is saved as ExtensionDensity_DublinED.png.
"""

from pathlib import Path

import duckdb
import geopandas as gpd
import matplotlib.pyplot as plt


#########################################
##### SETUP #####
#########################################

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
OUTPUTS = ROOT / "outputs"
OUTPUTS.mkdir(exist_ok=True)

PLANNING = DATA / "planning" / "planning_applications_national.gpkg"
EDS = DATA / "boundaries" / "eds_national.gpkg"
DATABASE = OUTPUTS / "extensions.duckdb"
DENSITY_FILE = OUTPUTS / "DublinED_ExtensionDensity.gpkg"
MAP_FILE = OUTPUTS / "ExtensionDensity_DublinED.png"

con = duckdb.connect(str(DATABASE))
con.install_extension("spatial")  # downloads it the first time only
con.load_extension("spatial")


#########################################
##### 1. IMPORT AND FILTER THE TABLES #####
#########################################

# Import the two national GeoPackages into DuckDB.
con.execute("CREATE OR REPLACE TABLE planning AS SELECT * FROM ST_Read('" + PLANNING.as_posix() + "')")
con.execute("CREATE OR REPLACE TABLE eds_national AS SELECT * FROM ST_Read('" + EDS.as_posix() + "')")

national_planning_count = con.execute("SELECT COUNT(*) FROM planning").fetchone()[0]
national_ed_count = con.execute("SELECT COUNT(*) FROM eds_national").fetchone()[0]

# Keep the EDs in Dublin City.
con.execute("""
CREATE OR REPLACE TABLE eds AS
SELECT * FROM eds_national
WHERE COUNTY_ENGLISH = 'DUBLIN CITY'
""")
ed_feature_count = con.execute("SELECT COUNT(*) FROM eds").fetchone()[0]

# Apply each planning filter and store the number of rows remaining.
con.execute("""
CREATE OR REPLACE TABLE dublin_applications AS
SELECT * FROM planning
WHERE PlanningAuthority = 'Dublin City Council'
""")
dublin_application_count = con.execute("SELECT COUNT(*) FROM dublin_applications").fetchone()[0]

con.execute("""
CREATE OR REPLACE TABLE dated_applications AS
SELECT * FROM dublin_applications
WHERE ReceivedDate >= '2019-01-01'
  AND ReceivedDate < '2026-01-01'
""")
dated_application_count = con.execute("SELECT COUNT(*) FROM dated_applications").fetchone()[0]

con.execute("""
CREATE OR REPLACE TABLE extensions AS
SELECT * FROM dated_applications
WHERE lower(DevelopmentDescription) LIKE '%extension%'
""")
extension_feature_count = con.execute("SELECT COUNT(*) FROM extensions").fetchone()[0]

# Report every filter count in two summary lines.
print("Planning rows: national =", national_planning_count,
      "; Dublin City Council =", dublin_application_count,
      "; 2019-2025 =", dated_application_count,
      "; extensions =", extension_feature_count)
print("ED rows: national =", national_ed_count,
      "; Dublin City =", ed_feature_count)


#########################################
##### 2. CHECK CRS AND REPROJECT #####
#########################################

# Store the geometry types for the layer summaries.
ed_geometry = con.execute("SELECT DISTINCT ST_GeometryType(SHAPE) FROM eds").fetchone()[0]
application_geometry = con.execute("SELECT DISTINCT ST_GeometryType(SHAPE) FROM extensions").fetchone()[0]
print("EDs:", ed_feature_count, "features; CRS = EPSG:2157; geometry =", ed_geometry)

# Reproject the points to Irish Transverse Mercator in SQL.
con.execute("""
CREATE OR REPLACE TABLE extensions_2157 AS
SELECT * REPLACE (
    ST_Transform(SHAPE, 'EPSG:3857', 'EPSG:2157') AS SHAPE
)
FROM extensions
""")
print("Applications:", extension_feature_count, "features; geometry =", application_geometry,
      "; CRS reprojected from EPSG:3857 to EPSG:2157")


#########################################
##### 3. SPATIAL JOIN AND COUNTS #####
#########################################

# Join each ED to the points within it and count the points.
con.execute("""
CREATE OR REPLACE TABLE extension_counts AS
SELECT
    d.ED_GUID,
    d.ED_ENGLISH,
    d.COUNTY_ENGLISH,
    d.SHAPE,
    COUNT(e.OBJECTID) AS extension_count
FROM eds d
LEFT JOIN extensions_2157 e ON ST_Within(e.SHAPE, d.SHAPE)
GROUP BY d.ED_GUID, d.ED_ENGLISH, d.COUNTY_ENGLISH, d.SHAPE
""")

# A left join in the opposite direction finds points outside every ED.
outside_count = con.execute("""
SELECT COUNT(*)
FROM extensions_2157 e
LEFT JOIN eds d ON ST_Within(e.SHAPE, d.SHAPE)
WHERE d.ED_GUID IS NULL
""").fetchone()[0]
joined_count = con.execute("SELECT SUM(extension_count) FROM extension_counts").fetchone()[0]
print("Applications joined to an ED =", joined_count,
      "; outside every ED =", outside_count)


#########################################
##### 4. CALCULATE DENSITY AND SAVE #####
#########################################

# EPSG:2157 uses metres, so divide square metres by one million.
con.execute("""
CREATE OR REPLACE TABLE extension_density AS
SELECT
    ED_GUID,
    ED_ENGLISH,
    COUNTY_ENGLISH,
    extension_count,
    ST_Area(SHAPE) / 1000000 AS area_km2,
    extension_count / (ST_Area(SHAPE) / 1000000) AS extension_density,
    SHAPE
FROM extension_counts
""")

# Remove an older output so that the GeoPackage can be replaced on a rerun.
DENSITY_FILE.unlink(missing_ok=True)
copy_sql = "COPY (SELECT * FROM extension_density) TO '"
copy_sql += DENSITY_FILE.as_posix()
copy_sql += "' WITH (FORMAT GDAL, DRIVER 'GPKG', LAYER_NAME 'DublinED_ExtensionDensity')"
con.execute(copy_sql)
print("Saved:", DENSITY_FILE)
con.close()


#########################################
##### 5. MAP EXTENSION DENSITY #####
#########################################

# Read the saved result and draw a simple choropleth map.
density = gpd.read_file(DENSITY_FILE)

fig, ax = plt.subplots(figsize=(9, 9))
density.plot(
    column="extension_density",
    cmap="YlOrRd",
    legend=True,
    edgecolor="grey",
    linewidth=0.3,
    ax=ax,
    legend_kwds={"label": "Extension applications per square kilometre"}
)
ax.set_title("Dublin City extension application density, 2019-2025")
ax.set_axis_off()
plt.savefig(MAP_FILE, dpi=200, bbox_inches="tight")
plt.close()
print("Saved:", MAP_FILE)
