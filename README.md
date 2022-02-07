# PCA on River networks

This repo contains the source code for PCA on river networks. Especially, it analyzes the Geum River TOC data. 


### Description

- Code
  - **source.R** and **iwanthue.R** include functions defined for analyses and visualization. 

  - **pca_river.R** is a main code for the data analysis.

- stpca_Rpackage
  - R package from "Flow-directed PCA (Gallacher et al., 2017)". 
  - It is used to implement PCA on spatiotemporal data with adjustment for spatial / temporal autocorrelation. 


### Code overview
The process described in **pca_river.R" is as follows.

- Data load / preprocessing
- Weight construction
- Application of flow-directed PCA 
- The proposed method
