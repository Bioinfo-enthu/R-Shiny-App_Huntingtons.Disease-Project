# R shiny app: Post-mortem Huntington’s Disease prefrontal cortex patient samples compared with neurologically healthy control samples

-> **App :**

R Shiny application with 4 tabs for analyzing mRNA-Seq expression profiling dataset of human post-mortem BA9 brain tissue. This app focuses on Huntington's Disease patient samples and included neurologically normal individuals samples as controls.

**1. Samples Tab:** Provided a comprehensive summary of samples. Fetched metadata from CSV files. Plotted bar graphs for different X-axis factors.

**2. Counts Tab:** Filtered normalized counts data based on gene variance and non-zero values. Implemented functionality to plot scatter plots, heatmaps, and PCA according to filter values.

**3. Differential Expression Tab:** Identified differentially expressed genes and presented them in a tabular format. Plotted a volcano plot to visually represent the differentially expressed genes.

**4. Correlation Expression Analysis Tab:** Established a correlation network between differentially expressed genes. Generated a clustered heatmap for visualization. Fetched correlation network metrics for the specified genes.

-> **Dataset Used :**

**mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington’s Disease and neurologically normal individuals**

**https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810**

This dataset profiled gene expression with RNASeq in post-mortem human dorsolateral prefrontal cortex from patients who died from Huntington’s Disease and age and sex matched neurologically healthy controls.

