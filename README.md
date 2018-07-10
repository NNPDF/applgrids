# applgrids
APPLgrids for NNPDF

## Project summary and aim

The aim of this repository is to store the "official" APPLgrids used by the NNPDF Collaboration.
Several grids stored here have been computed with the `external` repository, 
however, a smaller number of grids have been generated from experimental collaborations or th. colleagues.

Each folder matches the name of the dataset registered in the `apfelcomb` project database and 
reflects the official name used in NNPDF fits through the `buildmaster` implementation.

Each folder contains a specific README file with the summary information about the grid origin and usage.

### Testing code

A small testing code for generating predictions from APPLgrids + LHAPDF is
provided in `appltest`.

### Tag policy

This repository is usually tagged before the production of a new NNPDF release so tables have a history.

### Commit policy

The update of tables in this repository can be performed directly from developers. 
New tables must agree with the naming conventions used by other projects.
