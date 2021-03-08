function [log_likeli_subj, params0 ] = cal_log_likeli_optimal( params_subj, data_subj, n_trials )
% secretary problem. gain and loss can be based on a round or a trial??.
% round based first. ?? ??? ???...12/27/18
% nesxt step: power function. usually two or three parameter. 

sigma = params_subj;

numbers = 1:150;
params0( n_trials - 1, 1) = 1/2 * max(numbers);

for i = 1:(n_trials - 2)
    params0(-i+n_trials-1) = (((max(numbers)-params0(-i+n_trials))/max(numbers))...
        *(params0(-i+n_trials)+max(numbers))/2) + params0(-i+n_trials)/max(numbers)*params0(-i+n_trials);
end

log_likeli_subj = 0;

for i = 1:(n_trials - 1)
    
    doi = data_subj( data_subj(:,3) == i, : );
    
    if size(doi, 1) > 1
        data = doi( : , 4:5 );
        thresh_u = params0(i);
        log_likeli_temp = cal_likeli_SP( [thresh_u, sigma], data, 1 );
        log_likeli_subj = log_likeli_subj + log_likeli_temp;
    end
    
end

log_likeli_subj = -log_likeli_subj;


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


