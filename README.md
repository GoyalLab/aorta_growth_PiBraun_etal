Repository for the Paper:
Resolving the Design Principles that Control Postnatal Vascular Growth and Scaling

Overview

This repository contains the analysis code used in the study titled "Resolving the Design Principles that Control Postnatal Vascular Growth and Scaling." The code is organized by dataset type and analysis steps to facilitate easy access and understanding.

Directory Structure

Nuclei Segmentation
- Cellpose Segmentation
- Nuclei counting: nuclei/cellSize_nucNumber_downstream.ipynb
- nuclei statistics: cellShape/getNucleiStatistics/cellSize_getNucleiInformation.ipynb
  
Cell Shape Analysis
- AutomatedQuality Check of Segmentation: cellShape/cellSize/cellSize_automatedQualityCheckCellSizes.ipynb
- Major/minor axis and aspect ratio: cellSize/getCellSizes_downstream.ipynb
- Density of the whole stitched image: cellShape/density/cellSize_densityStitched.ipynb
- Density of cropped images: cellShape/density/cellSize_densityTiles.ipynb
- Density and Cell Number Along Axis: cellShape/nucleiAlongAxis/cellSize_densityNumberCellsAlongLengthAndCircumference.ipynb
- Adjust length, major axis, and cell area using microCT measurements: cellShape/adjust_microCT/cellSize_adjustLengthArea.ipynb
- Length and Circumference Estimation: cellShape/nucleiAlongAxis/cellSize_Length_CircumferenceEstimate.ipynb
  
Ki67 Active Analysis
- Ki67 Signal Extraction: ki67/extraction/ki67Active.ipynb
- Visualize and define thresholding for Ki67 active cells: ki67/extraction/ki67_visualizationAndThresholding.ipynb
- Applying Thresholding to Cells: ki67/thresholding/KI67_downstream_Applythresholding.ipynb
- 
Inducible Rainbow Analysis
- Watershed Signal Segmentation: inducibleRainbow/segmentation/waterShedSegmentation.py
- Match Inducible Rainbow with Nuclei Segmentation: inducibleRainbow/segmentation/inducibleRainbow_segmentation.ipynb
- Create Final CSV for Clone Sizes: inducibleRainbow/downstream/inducibleRainbow_createCountStatCSV.ipynb
- Regression for clonal analysis: inducibleRainbow/regression/indRainbow_anglesAndRegressionLineCellpose.ipynb
- K-NN analysis of inducible rainbow clones: inducibleRainbow/statistics/inducibleRainbow_kNNNegativePositiveControl.ipynb
- Diversity Statistics Calculation (e.g., Gini): inducibleRainbow/diversity_Statistics.ipynb
- 
Ki67 and Inducible Rainbow
- Segmentation of inducible raibow signal and Numpy File Creation: inducibleRainbowKI67/segmentation/inducedRainbowKI67_segmentationInducedRainbowKI67.ipynb
- Cluster Size Mapping and mapping Ki67 active cells: inducibleRainbowKI67/indRainbowKI67_mappingArea_ClusterSize.ipynb
- Distribution Distance Calculations and Simulations for number of ki67 active cells per cluster size: inducibleRainbowKI67/distance_simulation/indRainbowKI67_distance.ipynb
- Statistical Testing on Cluster Distances: inducibleRainbowKI67/cellCycle_statTest_confidence_pval_singletons.ipynb
- Image Creation for Active Clusters with matching colors for clustered cells and nuclei: inducibleRainbowKI67/cellCycle_filterCluster.ipynb
  
EdU Analysis
- Extract Active EdU Cells: EdU/cellCycle_extractActive.ipynb
- Filter Smooth Muscle Cells in Nuclei Segmentation: EdU/cellCycelEdUOnly_filerSmoothMuscleCells_backgroundNucleiSignal.ipynb
- Identify Touching EdU Active Cells: EdU/cellCycleEduOnly_touchingCells.ipynb
- Area Analysis of Active Singletons: EdU/pruning_cellAnalysis.ipynb
- L2-Distance Calculation among singletons to other clusters, between clusters, and between cells within a cluster: EdU/cellCycleEdu_distanceSingletons.ipynb
- Statistical Testing for Cluster Distribution: EdU/cellCycle_statTest_confidence_pval_singletons.ipynb
- Visualize Cells in EdU Clusters: EdU/cellCycle_filterCluster.ipynb
  
Fucci Mice Analysis
- Extraction for CDT1 & Geminin: See Ki67 extraction and thresholding methods provided above.

OSX-Cre: 
- Cell and Nuclei Segmentation and signal extraction: Follow the previously outlined processes for segmentation and signal extraction.

scRNA-seq analyis: 
singleCell_RNAseq/seurat5.R


Additional Information
•	Ensure that all required libraries are installed before running the scripts.
For any inquiries or issues, please reach out through the repository or contact the corresponding author.

