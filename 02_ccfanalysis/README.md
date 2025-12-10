# 02_ccfanalysis

Runs cortical connective field modelling based on V1 vertices fitting with data of V2 vertices. For exact model description, see SPOT/documentation/CorticalConnectiveFields.pdf.

## Directory structure

`run_model.py`: contains all process of cortical connective field modelling. The functions used in the code are saved in `ccf_model`.
Input: Benson's template, white matter surface, curv surface, functional image
Output: center voxel of connective field in V1 (v0), distance from the seeds along cortical surface (sigma), model fit (resuidual sum of squares) between V1 seed vertex and V2 vertex (rss), variance explainedr^2 (r)

 