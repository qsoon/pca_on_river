# PCA on River networks

This repo contains the source code for PCA on river networks. Especially, it analyzes the Geum River TOC data.

For the detail, please refer to our paper: [https://doi.org/10.1002/env.2753](https://doi.org/10.1002/env.2753).


## Description

- Code
  - `source.R` and `iwanthue.R` include functions defined for analysis and visualization. 

  - `pca_river.R` is a main code for the data analysis.

- Data
  - `KRF_3.0_Geumgang` contains shape files related to the Geum River such as catchment area shape, line shape, nodes. They can be downloaded from [Korean Reach File](http://water.nier.go.kr/web/gisKrf?pMENU_NO=89).
  - `ProcessedData` contains processed TOC data for Geum River. The original TOC data can be obtained from [Water Environment Information System](http://water.nier.go.kr/web/waterMeasure?pMENU_NO=2). 
  
- stpca_Rpackage
  - R package from "Flow-directed PCA (Gallacher et al., 2017)". 
  - It is used to implement PCA on spatiotemporal data with adjustment for spatial / temporal autocorrelation. 

- figures
  - Visualization results of the analysis. 

## Code overview
The process described in `pca_river.R` is as follows.

- Data load / preprocessing
- Weight construction
- Application of flow-directed PCA 
- The proposed method
