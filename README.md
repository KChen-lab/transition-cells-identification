# Identifying and characterizing transition cells in developmental processes from scRNA-seq data 
Decoding cellular state transitions is crucial for understanding complex biological processes in development and disease. While recent advancements in single-cell RNA sequencing offer insights into cellular trajectories, existing tools primarily study expressional rather than regulatory state shifts. We present a novel framework utilizing gene-pair correlations to detect transition cells from scRNA-seq data. Applying our approach to various contexts, including tissue regeneration, preinvasive lesions and humoral responses post-vaccination, reveals transition cells and their distinct gene expression profiles. Our study sheds light onto the underlying molecular mechanisms driving cellular state transitions enhancing our ability to identify therapeutic targets for disease interventions.

Our method models gene expressions using stochastic differential equations (SDEs), and identifies transition cells based on gene pair-wise correlation coefficients in 4 steps: <br />
Step 1: Data pre-processing <br />
Step 2: Selecting cells for Pearson’s correlation coefficients calculation <br />
Step 3: Top most variable genes selection <br />
Step 4: Transition index calculation <br />

![workflow](https://github.com/KChen-lab/transition-cells-identification/blob/main/images/workflow.png)

For more informations and mathematical details please refer to our manuscript
