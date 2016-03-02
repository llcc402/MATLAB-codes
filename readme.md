# The Smoothed Hierarchical Dirichlet Process(SHDP)

## Functions included:

1. data_generation: Generate synthetic data set with the following processure:

  1.1. Generate G0

  1.2. Generate G1,...,G20 smoothly

  1.3. Add noise to G1,...,G20. The noise is dirichlet distributed

  1.4. Generate data set according to G1,...,G20

2. discreternd: Generate discrete distributed data

3. dpDisrnd: Generate DP distributed distributions with a discrete base measure 

4. gem: Generate DP distributed distribution 

5. kl: Compute KL-divergence between two discrete distributions

6. sample_the_first: used in smoothSample

7. sample_the_rest: used in smoothSample

8. smoothSample: sampling a discrete distribution with base measure G0 and force it to have small symetric KL with a given distribution

9. symKL: Compute the symetric KL between two discrete distributions. Whenever the prob of any distribution is zero, we set it to be 1e-5

10. synthc_particle: Use particle filtering to sample distributions smoothly with synthetic data set

11. traditional_particle: Use particle filtering to sample distributions (DP distributed) with traditional methods.

12. trucBetarnd: Generate truncated Beta distributed data.

13. mykmeans: used in spectral clustering.

14. sym_cluster: the spectral clustering with symmetric normalized Laplacian.

## Scripts included:

1. change_bound: show how the sym KL for G_2,1 and G_2,2 changes with different bound B.

2. compare_kl: compare the sym KL with SHDP and HDP.

3. main: generate data with data_generate and sample G_1:20 from SHDP and HDP.

4. pami: the application to PAMI data set with SHDP.

5. pami_hdp: the application to PAMI data set with HDP.

## Data included:

1. the pami data set. The rows are corresponding to pami papers and the columns are corresponding to words.
