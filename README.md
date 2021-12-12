# APriD

This repository is the implementations of the paper: Adaptive Primal-Dual Stochastic Gradient Method for Expectation-constrained Convex Stochastic Programs.

The code is in Matlab. 
For QCQP problems, CVX is used to get the optimal solution. CVX is a Matlab-based modeling system for convex optimization. We download CVX from http://cvxr.com/cvx/download/. Particularly, we use the mosek solver in CVX, which needs a license. See http://cvxr.com/cvx/doc/mosek.html for details.

"run_all.m" includes how to run the experiments. Running it would produce all the results in the paper.

The folder "Neyman_Pearson" includes all code for the Neyman-Pearson Classification Problem. 
data_set includes the preprocessed data sets {'spambase', 'madelon'} ('gisette' is not uploaded beacuses it is larger than 25MB). 
iALM_NP.m is the deterministic method iALM which gives an optimal  solution.
NP_class_ApriD.m, NP_class_CSA.m, NP_class_MSA.m  are ApriD, CSA, MSA methods, respectively.
compare_algs_NP_class.m runs all the above methods and plot_compare_algs_NP_class.m plots the results for a given data set's filename in {'spambase', 'madelon', 'gisette'}.
select_parameters_NP_class.m runs ApriD method with different parameters and plot_select_parameters_NP_class.m plots the results for a given data set's filename in {'spambase', 'gisette'}.

The folder "QCQP_expect" includes all code for QCQP in Expectation Form.  
QCQP_expect_ApriD.m, QCQP_expect_CSA.m, QCQP_expect_MSA.m are ApriD, CSA, MSA methods, respectively.
compare_algs_QCQP_expect.m runs all the three methods and plot_compare_algs_QCQP_expect.m plots the results for a given n in {10, 200}.
select_parameters_QCQP_expect.m runs ApriD method with different parameters and plot_select_parameters_QCQP_expect.m plots the results for a given n in {10, 200}.

The folder "QCQP_scenar" includes all code for QCQP in finite-sum structured QCQP with many constraints.  
QCQP_scenario_ApriD.m, QCQP_scenario_CSA.m, QCQP_scenario_MSA.m, QCQP_scenario_PDSG_adp are ApriD, CSA, MSA, PDSG_adp methods, respectively.
compare_algs_QCQP_scenario.m runs all the three methods and plot_compare_algs_QCQP_scenario.m plots the results for a given n in {10, 200}.
 
After running "run_all.m", "running_time.txt" will be created and includes the computing time of each methods in the three examples.
