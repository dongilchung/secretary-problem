function log_likeli = cal_log_likeli_movingAverage( params, data, n_alternatives )

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

log_likeli = -sum( log( P_correct ));

if mu > 150 || mu < 20 
    log_likeli = log_likeli + 10000;
end


