function UBSDM_MCMC_exp( data_behavior, n_samples, EXP, repeat_number )

output_file = strcat('UBSDM_Exp',num2str(EXP),'_r',num2str(repeat_number),'.mat');


n_params = 4;
% const.prior = 4;
const.prior = 2; % this is to match the mode of inverse gamma distribution properly. 
thinning = 60;
prop_update = 200;
factor_update = 200;
lambda1 = 1; % 0 to 1
lambda2 = 3; % 0 >
prop_factor = 1e-7;

% data_behavior = [ subjNum, rnd, trialIndex, value, took_it ];
subjList = unique( data_behavior(:,1) );
n_subjs = length(subjList);

% generate curr_sample0
load('UBSDM_Exp4_r1.mat')
reform_sample = zeros( n_params*(n_subjs+2), 2000);
for i = 1:2000
    reform_sample(:,i) = model_to_sample( samples(1:end,i+2000), 4 );
end
% 

curr_sample = mean( reform_sample, 2 );
std_temp = [ 0.7, 0.7, 0.5, 0.03 ]';
curr_sample(end-3:end) = std_temp;

var_temp = var( reform_sample, [], 2 );
var_std = 0.25*std_temp.^2;
var_temp = [ var_temp(1:end-4); var_std ];

prop_cov = diag( var_temp );
    
%%%
% LB = [eps, eps, eps, -0.1]';
% UB = [150, 1, 50, 0.5]';
%%%

%
log_prior_curr = cal_log_prior(curr_sample, const, n_params, n_subjs);
log_like_curr = cal_log_like(curr_sample, data_behavior, n_params, n_subjs, subjList);
log_post_curr = log_like_curr + log_prior_curr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
samples = zeros(length(curr_sample), n_samples);

update_times = 1;
update_accept = 0;
sample_number = 1;
accept = 0;
reject = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a = GetSecs;

for i = 1:n_samples*thinning
    
    %a = GetSecs;curr, prop_factor, prop_cov)
    proposal_sample = sampling_proposal(curr_sample, prop_factor*prop_cov);
    
    log_prior_prop = cal_log_prior(proposal_sample, const, n_params, n_subjs);
    log_like_prop = cal_log_like(proposal_sample, data_behavior, n_params, n_subjs, subjList);
    log_post_prop = log_prior_prop + log_like_prop;
    
    acc = log_post_prop - log_post_curr;
    if isnan(log_post_prop)
        acc = -inf;
    end
    
    if log(rand) <= acc
        curr_sample = proposal_sample;
        log_post_curr = log_post_prop;
        %log_like_curr = log_like_prop;
        accept = accept + 1;
        update_accept = update_accept + 1;
    else
        %samples(i+1) = samples(i);
        reject = reject + 1;
    end
    
    if ~rem(i,thinning)
        samples(:,sample_number) = curr_sample;
        sample_number = sample_number + 1;
        disp(round([sample_number-1,100*accept/i]))
%         disp(log_post_curr);
%         disp(log_like_curr);
%         disp(curr_sample);
    end
        
    gamma1 = (1/(update_times)^(lambda1));
    gamma2 = lambda2*gamma1;
        
    if ~rem(i, thinning*prop_update)
        start_number = max(sample_number - 800, 1);
        vector_update = samples(:, start_number:(sample_number-1));
        prop_cov = cov(vector_update');
    end
    
    if ~rem(i,factor_update)
        log_prop_factor = log(prop_factor) + gamma2 * ((update_accept / factor_update) - 0.234);
        prop_factor = exp(log_prop_factor)
        disp(update_accept / factor_update)
        update_accept = 0;
    end
    
    %b = GetSecs;
    %disp([i, accept])
    %duration = b-a
    
end
disp(round([sample_number,100*accept/i]))
% b = GetSecs;
%disp([i, accept, b-a])
% duration = b-a

for i = 1:n_samples
    samples(:,i) = sample_to_model(samples(:,i), n_params);
end


save (output_file, 'samples');

% matlabpool close



function log_prior = cal_log_prior(sample, const, n_params, n_subjs)

sample_matrix = reshape( sample, n_params, length(sample)/n_params );

mu_vector = sample_matrix(:,end-1);
mu_vector = repmat( mu_vector, n_subjs, 1 );

std_prior = sample_matrix( :, end );
% std_prior = exp( std_prior );

std_vector = repmat( std_prior, n_subjs, 1 );

log_prior_all = log(normpdf( sample(1:n_params*n_subjs), mu_vector, std_vector ));
log_prior = sum(log_prior_all);

B = const.prior*[ 0.4, 0.4, 0.5, 0.02 ]';
A = 1;
% v = 1;

if sum( std_prior' <= [0.01, 0.01, 0.01, 0.01] ) > 0 
    prior_sigma = 0;
else
    % prior_sigma = ( 1 + (1/v)*( std_prior./A ).^2).^-(v+1)/2;
    prior_sigma = sum(log( inversegampdf_tune( std_prior, A, B ) ));
end

% log_prior = log_prior + log( prior_sigma ); Ooops! I applied log twice!! 
log_prior = log_prior + prior_sigma;


function sum_log_likeli = cal_log_like( curr_sample, data_all, n_params, n_subjs, subjList )

log_likeli_all = zeros( n_subjs, 1 );

params_vector = sample_to_model( curr_sample, n_params );
params_matrix = reshape( params_vector, n_params, n_subjs + 2 );

for subj = 1:n_subjs
    
    params = params_matrix( :, subj );
    data_subjs = data_all( data_all(:,1) == subjList(subj), :);
    
    log_likeli_all(subj) = cal_log_likeli_subj_wc( params, data_subjs);
    
end
sum_log_likeli = sum( log_likeli_all );


function [log_likeli_subj, threshold, exp_value] = cal_log_likeli_subj_wc( params_subj, data_subj )
% secretary problem. gain and loss can be based on a round or a trial??.
% round based first. ?? ??? ???...12/27/18
% nesxt step: power function. usually two or three parameter. 

reference = params_subj(1);
power_param = params_subj(2);
sigma = params_subj(3);
wc = params_subj(4); % waiting cost in utility

n_trials = 10;

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
    
    data = doi( : , 4:5 );
    thresh_u = threshold(i);
    
    log_likeli_temp = cal_likeli_SP( [thresh_u, sigma], data, 1 ); 
    log_likeli_subj = log_likeli_subj + log_likeli_temp;
   
end


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

function sample_all = sampling_proposal(curr_vector, prop_cov)

sample = mvnrnd( curr_vector', prop_cov );
sample_all = sample';


function model_params = sample_to_model( curr_sample, n_params )
% vector to vector
sample = reshape( curr_sample, n_params, length(curr_sample)/n_params );

ref_temp = sample( 1, 1:end-1 );
power_temp = sample( 2, 1:end-1 );
sigma_temp = sample( 3, 1:end-1 );

sample( 1, 1:end-1 ) = 150./(1+exp(-ref_temp));
sample( 2, 1:end-1 ) = 2./(1+exp(-power_temp));
sample( 3, 1:end-1 ) = exp( sigma_temp );
% sample( :, end ) = exp( sample( :, end ) );

model_params = sample(:);


function sample = model_to_sample( model_params, n_params )
% vector to vector

params_matrix = reshape( model_params, n_params, length(model_params)/n_params );

ref_temp = params_matrix( 1, 1:end-1 );
power_temp = params_matrix( 2, 1:end-1 );
sigma_temp = params_matrix( 3, 1:end-1 );

params_matrix( 1, 1:end-1 ) = log( ref_temp./(150 - ref_temp));
params_matrix( 2, 1:end-1 ) = log( power_temp./(2 - power_temp));
params_matrix( 3, 1:end-1 ) = log( sigma_temp );
% params_matrix( :, end ) = log( params_matrix( :, end ) );

sample = params_matrix(:);

