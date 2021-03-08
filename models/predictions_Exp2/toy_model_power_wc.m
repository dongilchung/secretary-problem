function [threshold, exp_value] = toy_model_power_wc( power_param, reference, wc, n_trials )
% secretary problem. gain and loss can be based on a round or a trial??.
% round based first. ?? ??? ???...12/27/18
% nesxt step: power function. usually two or three parameter. 

% reference = params_subj(1);
% power_param = params_subj(2);
% sigma = params_subj(3);
% wc = params_subj(4); % waiting cost in utility

% n_trials = 5;

min_v = 0;
max_v = 150;
step_size = 0.01;
x_vector = min_v:step_size:max_v;

exp_u_trial(n_trials) = cal_exp_util( min_v, max_v, power_param, reference, x_vector );
exp_u_cum(n_trials) = exp_u_trial(n_trials) - wc;
% exp_u_cum(n_trials) = exp_u_trial(n_trials);
exp_v_cum(n_trials) = 0.5*max_v;

for i = 1:(n_trials-1)
    
    index = n_trials - i;
    threshold(index) = trans_u_to_v( exp_u_cum(index + 1), reference, power_param, x_vector );
    exp_u_trial(index) = cal_exp_util( threshold(index), max_v, power_param, reference, x_vector );
    exp_u_cum(index) = ((max_v - threshold(index))/(max_v - min_v))*exp_u_trial(index) +...
        ((threshold(index)-min_v)/(max_v - min_v))*exp_u_cum(index + 1) - wc;
%     exp_u_cum(index) = ((max_v - threshold(index))/(max_v - min_v))*exp_u_trial(index) +...
%         ((threshold(index)-min_v)/(max_v - min_v))*exp_u_cum(index + 1);
    
    exp_v_trial(index) = 0.5*( threshold(index) + max_v );
    exp_v_cum(index) = ((max_v - threshold(index))/(max_v - min_v))*exp_v_trial(index) + ((threshold(index)-min_v)/(max_v - min_v))*exp_v_cum(index + 1);
    
end

exp_value = exp_v_cum(1);



function v_value = trans_u_to_v( u_value, exp_v_round, power_param, x_vector )

% u_value = u_value - wc;

u_x = cal_u_x( power_param, exp_v_round, x_vector );
[~, v_index] = min( abs(u_x - u_value));
v_value = x_vector(v_index);


function exp_util = cal_exp_util( min_v, max_v, power_param, reference, x_vector )

u_x = cal_u_x( power_param, reference, x_vector );

[~, min_v_index] = min( abs(x_vector - min_v));
[~, max_v_index] = min( abs(x_vector - max_v));

% exp_util = (1/(max_v-min_v))*(x_vector(2) - x_vector(1))*sum(u_x(min_v_index:max_v_index));
exp_util = mean(u_x(min_v_index:max_v_index));
% Can it be simply mean(u_x(min_v_index:max_v_index))?

