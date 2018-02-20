# spline-granger-causality
A procedure to reduce the parameters estimate in Granger-causal analysis through temporal smoothing.

Code to infer networks using a modified version of conditional Granger causality as described in <i> A procedure to increase power of Granger-causal analysis through temporal smoothing</i>.

Must have [MVGC Multivariate Granger Causality toolbox](http://users.sussex.ac.uk/~lionelb/MVGC/html/mvgchelp.html) to run as implemented in:

<i> Barnett, L., Seth, A.K., 2014. [The MVGC multivariate Granger causality toolbox: A new approach to Granger-causal inference](https://www.ncbi.nlm.nih.gov/pubmed/24200508). J. Neurosci. Methods 223, 50–68. doi:10.1016/j.jneumeth.2013.10.018. </i>.

Three example simulations are provided:

Example 1) main_sim_1N_high_freq.m
Example 2) main_sim_1N_low_freq.m
Example 3) main_sim_9N.m

Examples 1 and 2 fit models for one signal dominated by high or low frequencies, respectively.  Example 3 fits networks to a multivariate system of nine signals.

To execute an example, call the simualtion at the MATLAB command line:

```
>> main_sim_1N_high_freq
```

