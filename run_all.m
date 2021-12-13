clc
close all
clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run the results in Subsection 5.1-5.4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% creat the file to save the running time shown in Table 2
file_time = fopen('running_time.txt','w');
fprintf(file_time, ['Example' ' & ' 'data or size' ' & ' 'APriD' ' & ' 'MSA' ' & ' 'CSA' ' & '  'PDSG_adp'  ' \\\\\n']);
fclose(file_time);

%% run the code for Neyman-Pearson Classification Problem
addpath ./Neyman_Pearson
addpath ./Neyman_Pearson/data_set

fprintf('-----------------------------------------------------------------\n')
fprintf('run the code for Neyman-Pearson Classification  Problem\n')
fprintf('-----------------------------------------------------------------\n')
% run Neyman-Pearson Classification Problem with three different dataset
filename = 'spambase'; compare_algs_NP_class; clear; close all; fprintf('\n\n')
filename = 'madelon';  compare_algs_NP_class; clear; close all; fprintf('\n\n')
filename = 'gisette';  compare_algs_NP_class; clear; close all; fprintf('\n\n')

%% run the code for QCQP in Expectation Form
addpath ./QCQP_expect

fprintf('-----------------------------------------------------------------\n')
fprintf('run the code for QCQP in Expectation Form\n')
fprintf('-----------------------------------------------------------------\n')
n = 10;  compare_algs_QCQP_expect; clear; close all; fprintf('\n\n')
n = 200; compare_algs_QCQP_expect; clear; close all; fprintf('\n\n') % need a long time for this one

%% run the code for finite-sum structured QCQP with many constraints
addpath ./QCQP_scenar

fprintf('-----------------------------------------------------------------\n')
fprintf('run the code for finite-sum structured QCQP with many constraints\n')
fprintf('-----------------------------------------------------------------\n')
n = 10;  compare_algs_QCQP_scenario; clear; close all; fprintf('\n\n')
n = 200; compare_algs_QCQP_scenario; clear; close all; fprintf('\n\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% run the results in Subsection 5.5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run ApriD with different parameters for Neyman-Pearson Classification Problem
addpath ./Neyman_Pearson
addpath ./Neyman_Pearson/data_set

fprintf('-----------------------------------------------------------------\n')
fprintf('run ApriD with different parameters for Neyman-Pearson Classification Problem\n')
fprintf('-----------------------------------------------------------------\n') 
filename = 'spambase'; select_parameters_NP_class; clear; close all; fprintf('\n\n') 
filename = 'gisette';  select_parameters_NP_class; clear; close all; fprintf('\n\n')

%% run ApriD with different parameters for QCQP in Expectation Form
addpath ./QCQP_expect

fprintf('-----------------------------------------------------------------\n')
fprintf('run ApriD with different parameters for QCQP in Expectation Form\n')
fprintf('-----------------------------------------------------------------\n')
n = 10;  select_parameters_QCQP_expect; clear; close all; fprintf('\n\n')
n = 200; select_parameters_QCQP_expect; clear; close all; fprintf('\n\n') % need a long time for this one