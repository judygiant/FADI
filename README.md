# Reproducibility Files for FADI 

## R codes
The R codes folder contains the R scripts for simulation studies, and application of FADI to the 1000 Genomes data (estimation of principal eigenspace and inferential analysis under the degree-corrected mixed membership model). 

### Simulation Studies

The following scripts implement simulation experiments under various statistical models:
- `example_spiked_covariance.R`: Spiked covariance model  
- `example_GMM.R`: Gaussian mixture model (GMM)  
- `example_DCMM.R`: Degree-corrected mixed membership model (DCMM)  
- `example_missing_matrix.R`: Incomplete matrix inference model  

Each script accepts the following input arguments:
- `d`: Dimension of the data  
- `mc`: Index for independent Monte Carlo replicates  
- `rt`: Ratio of \( Lp/d \)

These arguments should be specified when invoking the script, for example via a shell script. Below is a sample SLURM batch script:

```bash
#!/bin/bash
#SBATCH --array=1-2
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -t 0-1:00
#SBATCH -p test
#SBATCH --mem=10G
#SBATCH -o %A_%a.out
#SBATCH -e %A_%a.err
#SBATCH --mail-type=NONE

module load gcc R  # Load R module
R CMD BATCH --quiet --no-restore --no-save "--args 1000 $SLURM_ARRAY_TASK_ID 1" example_GMM.R example_GMM_${SLURM_ARRAY_TASK_ID}.out
```
R scripts 1000g_estimation_layer_1.R and 1000g_estimation_layer_2.R contain the codes for applying FADI for estimating the principal eigenspace of the 1000 Genomes data. 1000g_estimation_layer_1.R implements Step 1 of FADI, with input parameters i-index for distributed data split, l-index for parallel sketching, and p-dimension of fast sketching. 1000g_estimation_layer_2.R implements step 2 of FADI, with input parameters l-index for parallel sketching, and p-dimension of fast sketching.

R script inference_1000g_SBM.R implements Step 1 and Step 2 of FADI for computing the top PCs of the undirected graph generated based on the 1000 Genomes data, with input parameter  l-index for parallel sketching. R script multiple_testing_1000g.Rmd performs multiple testing on inferring subject population of the 1000 Genomes data. The script multiple_testing_1000g.Rmd first implements Step 3 of FADI, by aggregating the parallel sketching results to output the FADI estimator of the top PCs. Then the script multiple_testing_1000g.Rmd performs the inferential procedure for membership testing using the FADI PC estimator, as detailed in Supplement D of the paper "Dimension Reduction for Large-Scale Federated Data: Statistical Rate and Asymptotic Inference".

## Data
The folder Data contains supplementary data for applying FADI to inferential analysis of the 1000 Genomes data. The file 1000g_sbm95.RData contains an undirected graph generated based on the 1000 Genomes data, used for the inferential application of FADI, and the file 1KG_TRACE_pca.txt contains the population information of the 1000 Genomes data subjects.

## Workflow
![FADI_workflow](FADI_workflow.png)
