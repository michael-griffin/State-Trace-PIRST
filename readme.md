## Repository for the data and analyses that were presented in:
### Benjamin, A.S., Griffin, M., & Douglas, J.V. & Dames, H. (2019). A nonparametric technique for analysis of state-trace functions


**Summary**  
State-trace analysis provides a direct way of evaluating a common question in studies of cognitive function: do one or two latent processes underlie performance on a particular task? In research on recognition memory for instance, items judged as seen before are often theorized as having an overall *familiarity* to them, in addition to a specific *recollection* of the time and circumstance it was last encountered. These theories compete with single process accounts: that a single latent *memory strength* measure is sufficient to account for the data collected. 

In this paper, we developed and benchmarked the PIRST technique, an algorithm based on isotonic regression that can be applied to state-trace plots and determine the amount of evidence present for >1 latent processes. The technique is benchmarked using simulated data, and is compared to an alternative approach, CMR analysis.

**Organization - Folders**

- **CMR analysis/** R code to run the 'coupled monotonic regression' analysis as described in *A statistical test of the equality of latent orders* by Kalish and colleagues.
- **CMR analysis/simulated data/** data simulated under varying conditions of noise, sparsity, and underlying shape, loaded by the CMR analysis
- **PIRST analysis/** Matlab code that generates the simulated data above, and runs the PIRST analysis
- **Simulation Results.../** Excel sheets with the test statistics of each analysis


