function [log_likeli_subj, threshold, exp_value] = cal_log_likeli_subj_wc( params_subj, data_subj, n_trials )
% secretary problem. gain and loss can be based on a round or a trial??.
% round based first. ?? ??? ???...12/27/18
% nesxt step: power function. usually two or three parameter. 

reference = params_subj(1);
power_param = params_subj(2);
sigma = params_subj(3);
wc = params_subj(4); % waiting cost in utility
% 
% if nargin == 2
%     n_trials = 5;
% end
% n_trials = max( data_


min_v = 0;
max_v = 150;
step_size = 0.01;
x_vector = min_v:step_size:max_v;

exp_u_trial(n_trials) = cal_exp_util( min_v, max_v, power_param, reference, x_vector );
exp_u_cum(n_trials) = exp_u_trial(n_trials) - wc;
exp_v_cum(n_trials) = 0.5*max_v;

for i = 1:(n_trials-1)
    
    index = n_trials - i;
    threshold(index) = trans_u_to_v( exp_u_cum(index + 1), reference, power_param, x_vector );
    exp_u_trial(index) = cal_exp_util( threshold(index), max_v, power_param, reference, x_vector );
    exp_u_cum(index) = ((max_v - threshold(index))/(max_v - min_v))*exp_u_trial(index) +...
        ((threshold(index)-min_v)/(max_v - min_v))*exp_u_cum(index + 1) - wc;
    
    exp_v_trial(index) = 0.5*( threshold(index) + max_v );
    exp_v_cum(index) = ((max_v - threshold(index))/(max_v - min_v))*exp_v_trial(index) + ((threshold(index)-min_v)/(max_v - min_v))*exp_v_cum(index + 1);
    
end

exp_value = exp_v_cum(1);

log_likeli_subj = 0;

for i = 1:(n_trials-1)
    
    doi = data_subj( data_subj(:,3) == i, : );
    
    %     util = cal_u_x( power_param, reference, doi(:,4) );
    %     thresh_u = cal_u_x( power_param, reference, threshold(i) );
    %     data = [ util', doi(:,5) ];
    
    if size(doi, 1) > 1
        data = doi( : , 4:5 );
        thresh_u = threshold(i);
        log_likeli_temp = cal_likeli_SP( [thresh_u, sigma], data, 1 );
        log_likeli_subj = log_likeli_subj + log_likeli_temp;
    end
    
end

log_likeli_subj = -log_likeli_subj;

% if wc < 0 || power_param > 1 || reference < 0 || reference > 150 
% if power_param > 1 || reference < 0 || reference > 150 
%     log_likeli_subj = log_likeli_subj + 1000;
% end




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


function log_likeli = cal_likeli_SP( params, data, n_alternatives )

durations = data(:,1);

mu = params(1);
sigma = params(2);

correct_or_not = data(:,2);

P_cdf = 0.5*erfc( -(durations-mu)./(sqrt(2)*sigma));% normcdf( durations, mu, sigma );
if n_alternatives == 2
    P_correct = 0.5 + 0.5*P_cdf;
elseif n_alternatives == 4
    P_correct = 0.25 + 0.75*P_cdf;
elseif n_alternatives == 1
    P_correct = P_cdf;
end
    
P_correct( P_correct == 1 ) = 1-eps;
P_correct( P_correct == 0 ) = eps;

P_correct( correct_or_not == 0 ) = 1 - P_correct( correct_or_not == 0 );

log_likeli = sum( log( P_correct ));


