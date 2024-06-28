<!-- title: Identifying Spatially Variable Genes  -->
# Identifying Spatially Variable Genes using Big-Small Patch (BSP) algorithm and sparse matrix implementation as scBSP

## BSP Algorithm

Big-small patch (BSP) is a granularity-guided, data-driven, and parameter-free model for identifying spatial variable genes in 2D and 3D high-throughput spatial transcriptomics data.

![BSP](https://github.com/juexinwang/BSP/blob/main/flowchart.png)

Original implementation of BSP in python is available at https://github.com/juexinwang/BSP


## scBSP - A Fast Tool for Single-Cell Spatially Variable Genes Identifications on Large-Scale Spatially Resolved Transcriptomics Data


This package utilizes a granularity-based dimension-agnostic tool, single-cell big-small patch (scBSP), implementing **sparse matrix** operation and KD-tree/balltree method for distance calculation, for the identification of spatially variable genes on
large-scale data. A corresponding Python library is available at [https://pypi.org/project/scbsp](https://pypi.org/project/scbsp/). The R source code is available at https://github.com/CastleLi/scBSP
. The Python source code is available at https://github.com/YQ-Wang/scBSP


# System Requirement
Tested on MacOS Sonoma version 14.4.1 with R version 4.3.1, and on Windows 11 12th Gen Intel(R) i7-1265U with R version 4.2.1.

# Installation
This package can be installed on R CRAN
```
install.packages("scBSP")
```
For installation errors with 'sparseMatrixStats', install the bioconductor package from [here](https://www.bioconductor.org/packages/release/bioc/html/sparseMatrixStats.html) before installing scBSP.

# Example 1: Breast Cancer Data
In your working directory, create a 'data' subdirectory to download the tutorial data files into. Load the scBSP package and Breast cancer data set, which can be downloaded [here](https://github.com/juexinwang/Tutorial_DahShu2024/blob/master/data/Layer2_BC_Count.rds). 

```
> library(scBSP)
> load("data/Layer2_BC_Count.rds")
```
     
View the Breat Cancer expression count matrix 'rawcount', where each row denotes a gene and each column represents a cell/spot.

```
> rawcount[1:5,1:5]

    17.907x4.967   18.965x5.003   18.954x5.995    17.846x5.993 20.016x6.019
GAPDH           1             7              5               1            2
USP4            1             0              0               0            0
MAPKAPK2        1             1              0               0            1
CPEB1           0             0              0               0            0
LANCL2          0             0              0               0            0
```

## Prepare the Input: Coordinates and Gene Expression

Extract the coordinates from the raw data
```
> info <- cbind.data.frame(
    x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)), 
    y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2))
    )
> rownames(info) <- colnames(rawcount)
> Coords <- info[,1:2]
```

View the coordinates of x and y on 2-D space

```
> head(Coords)

                  x     y
17.907x4.967 17.907 4.967
18.965x5.003 18.965 5.003
18.954x5.995 18.954 5.995
17.846x5.993 17.846 5.993
20.016x6.019 20.016 6.019
20.889x6.956 20.889 6.956
```

Check the coordinates dimensions (Should be M x D, representing D-dimensional coordinates for M spots)

```
> dim(Coords)
[1] 251   2
```

## Gene Expression Preprocessing
Exclude low expressed genes from the expression matrix
```
> Filtered_ExpMat <- SpFilter(rawcount)
```

For the gene expression data, we need to explicitly convert the data format to a sparse Matrix
```
> Filtered_ExpMat_sparse <- Matrix::Matrix(Filtered_ExpMat, sparse = TRUE)
```

Check the dimensions of the filtered input data (Should be N x M, sparse matrix with N genes and M spots)
```
> dim(Filtered_ExpMat_sparse)
[1] 11769   251
```

## Computing P-values
The inputs are the coordinates and the expresson matrix, scBSP calculates the p-values
```
> P_values <- scBSP(Coords, Filtered_ExpMat_sparse)
```

## Output the Final Results
```
> head(P_values)

  GeneNames     P_values
1     GAPDH 2.704573e-09
2      USP4 2.452562e-01
3  MAPKAPK2 4.177737e-03
4     CPEB1 7.656481e-01
5    LANCL2 7.272278e-01
6      MCL1 9.308317e-06
```

Display significant SVGs by sorting and filtering for P_values < 0.05 
```
> results <- P_values[P_values$P_values<0.05,]
> results <- results[order(results$P_values),]
> head(results)

     GeneNames     P_values
326     SPINT2 0.000000e+00
261      POSTN 6.439294e-15
2013    COL1A1 6.550316e-15
1347       B2M 7.660539e-15
2370       FN1 1.521006e-14
1072    EEF1A1 1.356693e-13

> dim(results)
[1] 1765    2
```
We can see there is a significant difference in the SVGs found by scBSP which have a very high P_value (PCLO) & genes which were found to have very low P_values (POSTN)

| PCLO  | POSTN |
| ------------- | ------------- |
| ![PCLO](plots/PCLO.png)  | ![POSTN](plots/POSTN.png) |


# Example 2: HDST data of the Mouse Hippocampus 
HDST provides a subcellular resolution spatial trianscriptomics data, which contains 181,367 spots and 19,950 genes with a density of ~0.0003. The data can be downloaded at [here](https://github.com/juexinwang/Tutorial_DahShu2024/blob/master/data/CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds)

```
> library(scBSP)
> load("data/CN24_D1_unmodgtf_filtered_red_ut_HDST_final_clean.rds")
```

View the expression count matrix which is already loaded in as a sparse matrix, each row denotes a gene and each column represents a cell/spot.

```
> head(sp_count)
Loading required package: Matrix
6 x 181367 sparse Matrix of class "dgCMatrix"
  [[ suppressing 34 column names ‘1000x100’, ‘1000x103’, ‘1000x113’ ... ]]
                                                                           
Rcn2    . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Mycbp2  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
mt-Rnr2 . . . . . . . . . . 1 . 1 2 . . 1 . 1 1 . . 2 2 2 6 4 . 6 2 . 1 1 4
Mprip   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Mroh1   . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
Zfp560  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
              
Rcn2    ......
Mycbp2  ......
mt-Rnr2 ......
Mprip   ......
Mroh1   ......
Zfp560  ......

 .....suppressing 181333 columns in show(); maybe adjust options(max.print=, width=)
 ..............................
```

## Prepare the Input: Coordinates and Gene Expression

Extract the coordinates
```
> info <- cbind.data.frame(
    x=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",1)),
    y=as.numeric(sapply(strsplit(colnames(sp_count),split="x"),"[",2))
    )
> rownames(info) <- colnames(sp_count)
> Coords <- as.matrix(info)
```

View the coordinates of x and y on 2-D space

```
> head(Coords)

            x   y
1000x100 1000 100
1000x103 1000 103
1000x113 1000 113
1000x114 1000 114
1000x116 1000 116
1000x142 1000 142
```

Check the coordinates dimensions (Should be M x D, representing D-dimensional coordinates for M spots)

```
> dim(Coords)
[1] 181367      2
```

## Gene Expression Preprocessing
Remove mitochondrial genes from the expression matrix
```
> mt_idx <- grep("mt-",rownames(sp_count))
> if(length(mt_idx)!=0){
    sp_count <- sp_count[-mt_idx,]
}
```

Excluding low expressed genes from the expression matrix
```
> Filtered_ExpMat <- SpFilter(sp_count)
```

Check the dimensions of the filtered input data (Should be N x M, sparse matrix with N genes and M spots)
```
> dim(Filtered_ExpMat)
[1]  13209 181367
```

## Computing P-values
Calculate Spatially Variable Genes with the filtered sparse matrix and coordinates using scBSP
```
> P_values <- scBSP(Coords, Filtered_ExpMat)
```

## Recording scBSP's Computation Time
Install peakRAM package to record computational resources
```
> install.packages("peakRAM")
```
Alternate installation method

```
> install.packages("remotes")
> remotes::install_github("tpq/peakRAM")
```

Record the computational resources used when running scBSP (It takes seconds on a laptop)
```
> library(peakRAM)
> peakRAM(P_values <- scBSP(Coords, Filtered_ExpMat))

Normalizing the expression matrix
Normalizing the coords matrix
Calculating p-values
                            Function_Call
1 P_values<-scBSP(Coords,Filtered_ExpMat)
  Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
1             4.64                0.1             771.4
```

Calcuate the density of the data
```
> Matrix::nnzero(sp_count)
[1] 1217700

> dim(sp_count)
[1]  19913 181367

> dimNum<-dim(sp_count)
> dimNum[1]*dimNum[2]
[1] NA
Warning message:
In dimNum[1] * dimNum[2] : NAs produced by integer overflow

> Matrix::nnzero(sp_count)/dimNum[1]/dimNum[2]
[1] 0.0003371672
```
The density of the data is ~0.0003, and by using this sparse matrix format, can accelerate the computational processing by a lot.

## Output the Final Results

Display significant SVGs by sorting and filtering for P_values < 0.05 
```
> results <- P_values[P_values$P_values<0.05,]
> results <- results[order(results$P_values),]
> head(results)

    GeneNames     P_values
24    Gm42418 0.000000e+00
70      Cmss1 0.000000e+00
74     Camk1d 3.330669e-16
49       Gphn 7.771561e-16
64 CT010467.1 1.975087e-13
43       Cdk8 1.445843e-12

> dim(results)
[1] 1829    2
```

# Example 3: SlideSeq V2 data on Mouse Olfactory Bulb

Use SlideSeq V2 data with higher density.

## Install data using SeuratData
```
> install.packages("remotes")
> remotes::install_github("satijalab/seurat-data")

> library(Seurat)
> library(SeuratData)
> library(SeuratObject)
> library(peakRAM)
```

Check the available data
```
> AvailableData()
```

Install Slide-seq v2 dataset of mouse hippocampus
```
> InstallData("ssHippo")
```

## Loading data and Preprocessing
Once installed, load the data into your enviorment and extract the filtered expression matrix
```
> slide.seq <- LoadData("ssHippo")
> data_extracted <- scBSP::LoadSpatial(slide.seq)
> ExpMatrix_Filtered <- scBSP::SpFilter(data_extracted$ExpMatrix, Threshold = 1)
```

Check the dimensions of the SlideSeq coordinates and filtered gene expression matrix
```
> dim(data_extracted$Coords)
[1] 53173     2

> dim(ExpMatrix_Filtered)
[1] 23243 53173
```

Record the computational resources used when running scBSP
```
> peakRAM({ P_values <- scBSP::scBSP(data_extracted$Coords,ExpMatrix_Filtered)})

Normalizing the expression matrix
Normalizing the coords matrix
Calculating p-values
                                                       Function_Call
1 {P_values<-scBSP::scBSP(data_extracted$Coords,ExpMatrix_Filtered)}
  Elapsed_Time_sec Total_RAM_Used_MiB Peak_RAM_Used_MiB
1             22.4                0.1            3485.7
```
We can see that it takes more time to run compared to the HDST data since they have different densities.

Calcuate the density of the data
```
> Matrix::nnzero(data_extracted$ExpMatrix)
[1] 22361774

> dim(data_extracted$ExpMatrix)
[1] 23264 53173

> dimNum<-dim(data_extracted$ExpMatrix)
> dimNum[1]*dimNum[2]
[1] 1237016672

> Matrix::nnzero(data_extracted$ExpMatrix)/(dimNum[1]*dimNum[2])
[1] 0.01807718
```
The density of the SlideSeq V2 data is ~0.02. The sparse matrix format cannot accelerate the computational time much as the HDST data which has a density of ~0.0003.

# Example 4: 3D Spatial data using Simulation

Since scBSP is dimension agnostic, we are able to use data of any dimension. Load the scBSP package and 3D simulated data file, which can be downloaded [here](https://github.com/juexinwang/Tutorial_DahShu2024/blob/master/data/Pattern_1.csv). 

## Install data using SeuratData
```
> library(scBSP)
> simdata <- read.csv("data/Pattern_1.csv")
```


## Prepare the Input: Coordinates and Gene Expression

Extract the coordinates from the simulated data file
```
> coords <- simdata[c(1,2,3)]
> colnames(coords) <- c("x", "y", "z")
```

View the coordinates of x, y, z in the 3-D space

```
> head(coords)

          x         y z
1 12.714909 13.352257 1
2  2.823547  4.759709 1
3 12.577456  3.674638 1
4 13.225454  2.316355 1
5  9.251169  4.510353 1
6  8.807597  7.141753 1
```

Check the coordinates dimensions (Should be M x D, representing D-dimensional coordinates for M spots)

```
> dim(coords)
[1] 2229   3
```

Extract that simulated gene counts and convert the format to a sparse matrix
```
> exp_matrix <- simdata[c(-1,- 2,-3)]
> exp_matrix <- Matrix::Matrix(t(exp_matrix), sparse = TRUE)
```

Check the dimensions of the filtered input data (Should be N x M, sparse matrix with N genes and M spots)
```
> dim(exp_matrix)
[1]   10 2229
```

## Computing P-values
Using the 3D coordatinates and expression matrix, scBSP calculates the p-values
```
> P_values <- scBSP::scBSP(coords,exp_matrix)
```

## Output the Final Results
```
> head(P_values)

  GeneNames   P_values
1     SVG_1 0.06804437
2     SVG_2 0.71548891
3     SVG_3 0.87505115
4     SVG_4 0.25470892
5     SVG_5 0.40861395
6     SVG_6 0.91673636
```

# Extra: Plotting Gene Expression in 2D Spatial Domain in Example 1
Additionally, to check the gene expression across the sample spots/cells, you can plot the gene expression with the spatial coordinates to visualize the distribution. Create a subdirectory in your working directory names 'plots' to save the figures.
(Note: Example 2 & 3's data is very large and may take up alot of computational resources to plot all the cells/spots. Please only use the below code for Example 1)

```
install.packages("ggplot2")
library(ggplot2)

saveGenePlotLog <- function(counts_matrix, Coords, geneName) {
  counts <- data.frame(t(data.frame(counts_matrix)))
  exp <- counts[,geneName]
  expo <- log(exp + 1)
  df1 <- data.frame(
    x = Coords$x,
    y = Coords$y,
    z = expo
  )
  p <- ggplot(df1, aes(x = x, y = y, color = z)) +
    geom_point(size = 4) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(color = geneName) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  ggsave(paste("plots/", geneName, ".png", sep = ""), plot = p, width = 8, height = 6, dpi = 300)
}
```
Once the function is defined, pass in the sample expression matrix, coordinates, and the name of the gene you want to create a plot for. In the example below, we use the expression matrix and coordinates from Example 1, and plot the SVG which was found to have the lowest P_value.
```
saveGenePlotLog(Filtered_ExpMat, Coords, "SPINT2")
```
| SPINT2  | 
| ------------- |
| ![SPINT2](plots/SPINT2.png)  | 

# Cite
1. Wang, J., Li, J., Kramer, S.T. et al. Dimension-agnostic and granularity-based spatially variable gene identification using BSP. Nat Commun 14, 7367 (2023). https://doi.org/10.1038/s41467-023-43256-5
2. Li, J., Wang, Y., Raina, M.A. et al. scBSP: A fast and accurate tool for identifying spatially variable genes from spatial transcriptomic data. bioRxiv 2024.05.06.592851; doi: https://doi.org/10.1101/2024.05.06.592851

# References:
1. https://github.com/juexinwang/BSP/
2. https://github.com/CastleLi/scBSP/
3. Stahl, P. L. et al. Visualization and analysis of gene expression in tissue sections by spatial transcriptomics. Science 353, 78-82, 2016
4. https://github.com/mssanjavickovic/3dst
5. https://xzhoulab.github.io/SPARK/02_SPARK_Example/

