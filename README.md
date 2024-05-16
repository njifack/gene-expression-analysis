# Gene Expression Analysis Project

This project involves the analysis of gene expression data to identify patterns and clusters related to different types of cancer. The analysis includes hierarchical clustering, k-means clustering, and k-medoids clustering methods to explore the relationships between gene expression profiles and cancer types.

### Overview

The project aims to:

- Identify gene expression patterns associated with different types of cancer.
- Cluster genes based on their expression profiles to uncover potential biomarkers.
- Evaluate the effectiveness of different clustering algorithms in grouping similar genes and cancer types.

### Data

The gene expression data used in this project is sourced from the Golub et al. (1999) study, which includes gene expression profiles for patients with acute lymphoblastic leukemia (ALL) and acute myeloid leukemia (AML). The dataset comprises gene expression levels for oncogenes and antigens.

### Analysis Methods

#### Hierarchical Clustering

Hierarchical clustering is performed using single and complete linkage methods to group genes based on their expression profiles. The resulting dendrograms visualize the hierarchical relationships between genes.

#### K-means Clustering

K-means clustering is applied to identify clusters of genes with similar expression patterns. The optimal number of clusters is determined by evaluating the sum of squared errors (SSE) for different values of k.

#### K-medoids Clustering

K-medoids clustering, also known as partitioning around medoids (PAM), is utilized to cluster genes based on dissimilarity measures. The clustering results are compared with those obtained from k-means clustering to assess the effectiveness of each method.

### Results and Interpretation

The clustering analyses reveal distinct groups of genes associated with different types of cancer. By examining the clusters and their gene composition, potential biomarkers for specific cancer types can be identified. Additionally, the comparison between hierarchical, k-means, and k-medoids clustering methods provides insights into their respective strengths and limitations in clustering gene expression data.

## How to Use

To replicate the analysis:

1. Load the gene expression data.
2. Perform hierarchical clustering using single and complete linkage methods.
3. Apply k-means clustering with different values of k and evaluate SSE.
4. Conduct k-medoids clustering with appropriate dissimilarity measures.
5. Compare the clustering results and interpret the findings in the context of cancer biology.

## Dependencies

- R programming language
- Required R packages: `multtest`, `cluster`, `FNN`

## Credits

- Golub et al. (1999) for providing the gene expression dataset.
- R Core Team and contributors for developing the R programming language and associated packages.

