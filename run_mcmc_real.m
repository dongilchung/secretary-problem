clear
clc
%
load('data_all_exp2_10.mat')

n_subj = size( data_all, 1 );
n_rounds = size( data_all, 2 );
%
data_behavior = format_data( data_all );
clear data_all

EXP = 4;
repeat_number = 1;
NSamples = 4000;
UBSDM_MCMC_exp( data_behavior, NSamples, EXP, repeat_number );
% UBSDM_MCMC_exp_p3( data_behavior, NSamples, EXP, repeat_number );
