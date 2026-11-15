# Sprint Log

## Sprint 0 (Set up)
Date： 2026-10-24

What was done:
Created the GitHub organization and the initial repository.
Set up the basic folder structure: data/, scripts/, and outputs/.
Downloaded the national planning applications dataset and the Electoral Divisions boundaries dataset.
Drafted the README.md  explaining the project overview, data sources, licences, and CRS.
Added the three script templates (Q1_ExtensionDensity.py, Q2_GlobalMoran.py, Q3_LISA_GiStar.py) into the scripts.
Discussed and assigned initial roles for the upcoming sprints.

What is next (Sprint 1):
Begin working on Question 1: Setting up the DuckDB spatial database, importing the data, and calculating extension density by ED.

## Sprint 1
Date:2026-11-03

## Coordinator Responsibilities
1.Review the progress and outputs from Sprint 1.
2.Confirm that the question 1 extension density Geopackage is available.
3.Coordinate the completion of Question 2.
4.Review team contributions and GitHub commits during today's coordination session.
5.Update the sprint log with progress, issues, and next steps.

## Issues and Dependencies

1.Question 1 script is available, but its generated outputs are not yet present in the GitHub repository.
2.Confirm availability of DublinED_ExtensionDensity.gpkg before running Question 2.

## Sprint 1 Review
2026-11-05
The question 1 was checked and runs correctly. The extensions.duckdb and the DublinED_ExtensionDensity.gpkg outputs are created and stored in the folder. The map is successfully created and saved in the output folder as well.
A line was added to the gitignore to avoid the push of the duckdb to Github which is too large.

### Sprint 2 
Date: 2026-11-10

Checked Q2 outputs and pulled the latest Q1 data for the analysis.
Computed Getis-Ord Gi* with 9999 permutations and applied FDR correction to identify hotspots and coldspots.
Merged the LISA script and Gi* code into `Q3_LISA_GiStar.py`, successfully resolving a Git merge conflict.
Generated the Gi* map and a cross-tabulation to compare LISA clusters with Gi* data.
