% predict the treshold of the 1st op of 2o experiment.

clear
clc

[~, exp_value5] = toy_model_power( 1, 100, 5 );
[~, exp_value2] = toy_model_power( 1, 100, 2 );
[~, exp_value10] = toy_model_power( 1, 100, 10 );

% exp_diff = exp_value5 - exp_value2;


% load('best_40_Exp2_p3.mat')
load('best_40_Exp2.mat')
% best_all = best_40_Exp2_p3;
best_all = best_40_Exp2;

n_repeat = 5000;
thresh_mean_all = zeros(1,n_repeat);

n_trials = 2;

exp_diff = -30:2:30;

for exp_diff_index = 1:length( exp_diff )
    
    exp_diff_current = exp_diff( exp_diff_index )
    
    for r = 1:n_repeat
        
        r;
        in_list = randi( [1, length(best_all)], 1, 20 );
        
        for i = 1:length(in_list)
            
            index = in_list(i);
            
            exp_v_round = best_all(index,1) + exp_diff_current; % expected difference from optimal models.
            loss_gain = best_all(index,2);
            
            
            [threshold, exp_value] = toy_model_power( loss_gain, exp_v_round, n_trials );
            
%             exp_v_diff = (exp_v_round - exp_value)^2;
%             
%             exp_value_all(i) = exp_value;
            threshold_all(i) = threshold(n_trials-1);
            
        end
        
        thresh_mean_all(r) = mean(threshold_all);
        
    end
    
    sorted  = sort( thresh_mean_all );
    min_boundary(exp_diff_index) = sorted(round(0.025*n_repeat));
    max_boundary(exp_diff_index) = sorted(round(0.975*n_repeat));
    
end
