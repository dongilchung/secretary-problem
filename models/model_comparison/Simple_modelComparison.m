clearvars; close all;

clear

index_4p = 1;
index_3p = 1;
index_linear = 1;
index_Olinear = 1;
index_const = 1;
index_opt = 1;
index_thr = 1;
index_OC = 1;

for Exp = 1

switch Exp
    
    case 1
        
        load('data_all_exp1.mat')
        subj_num = size( data_all, 1 );
        data_behavior = format_data( data_all );
        clear data_all
        load('best_40_Exp1.mat')
        best_params_all_p4 = best_40_Exp1;
        load('best_40_Exp1_p3.mat')
        best_params_all_p3 = best_40_Exp1_p3;
        load('thr_exp1.mat')
        thr = thr;
                
    case 2
        
        load('data_all_exp2_5.mat')
        subj_num = size( data_all, 1 );
        data_behavior = format_data( data_all );
        clear data_all
        load('best_40_Exp2_5.mat')
        best_params_all_p4 = best_40_Exp2;
        load('best_40_Exp2_5_p3.mat')
        best_params_all_p3 = best_40_Exp2_p3;
        load('thr_exp2_5.mat')
        thr = thr5;
        
    case 3
        
        load('data_all_exp3.mat')
        subj_num = size( data_all, 1 );
        data_behavior = format_data( data_all );
        clear data_all
        load('best_40_Exp3.mat')
        best_params_all_p4 = best_40_Exp3_12;
        load('best_40_Exp3_p3.mat')
        best_params_all_p3 = best_40_Exp3_p3;
        load('thr_exp3.mat')
        thr = thr;
        
    case 4
        
        load('data_all_exp2_10.mat')
        subj_num = size( data_all, 1 );
        data_behavior = format_data( data_all );
        clear data_all
        load('best_40_Exp2_10.mat')
        best_params_all_p4 = best_40_Exp4_2;
        load('best_40_Exp2_10_p3.mat')
        best_params_all_p3 = best_40_Exp4_p3;
        load('thr_exp2_10.mat')
        thr = thr;
        
end

% n_opp

n_trials = max( data_behavior(:,3));
n_data_points = size( data_behavior, 1 );

% cal likelihood of 4 params model

best_params_all = best_params_all_p4;
log_likeli_4p = 0;

for i = 1:subj_num
    
    data_subj = data_behavior( data_behavior(:,1) == i, : );
    best_params = best_params_all(i,:);
    
    [neg_log_likeli_subj, threshold, exp_value] = cal_log_likeli_subj_wc( best_params, data_subj, n_trials );
    log_likeli_subj = -neg_log_likeli_subj;
    
    % thresh_all(i,:) = threshold;
    log_likeli_all_4p(i) = log_likeli_subj;
    log_likeli_4p = log_likeli_4p + log_likeli_subj;
    
    n_data_points_subj = size( data_subj, 1 );
    
    AIC_4p_subj(index_4p) = 2*(4) -2*log_likeli_subj;
    BIC_4p_subj(index_4p) = log(n_data_points_subj)*4 -2*log_likeli_subj;
    index_4p = index_4p + 1;
end

AIC_4p = 2*(subj_num*4) -2*log_likeli_4p;
BIC_4p = log(n_data_points)*subj_num*4 -2*log_likeli_4p;


% cal likelihood of 3 params model

best_params_all = best_params_all_p3;
log_likeli_3p = 0;


for i = 1:subj_num
    
    data_subj = data_behavior( data_behavior(:,1) == i, : );
    best_params = best_params_all(i,:);
    
    [neg_log_likeli_subj, threshold, exp_value] = cal_log_likeli_subj( best_params, data_subj, n_trials );
    log_likeli_subj = -neg_log_likeli_subj;
    
    % thresh_all_p3(i,:) = threshold;
    log_likeli_all_3p(i) = log_likeli_subj;
    log_likeli_3p = log_likeli_3p + log_likeli_subj;
    
    n_data_points_subj = size( data_subj, 1 );
    AIC_3p_subj(index_3p) = 2*(3) -2*log_likeli_subj;
    BIC_3p_subj(index_3p) = log(n_data_points_subj)*3 -2*log_likeli_subj;
    index_3p = index_3p +1;
    
end

AIC_3p = 2*(subj_num*3) -2*log_likeli_3p;
BIC_3p = log(n_data_points)*subj_num*3 -2*log_likeli_3p;
% cal likelihood of optimal model

log_likeli_opt = 0;
params_subj0 = 20;

for i = 1:subj_num
    
    data_subj = data_behavior( data_behavior(:,1) == i, : );
    
    best_params = fminsearch( @(params)...
        cal_log_likeli_optimal( params, data_subj, n_trials ), params_subj0 );
    
    % best_params = [120, 0.2077, 0.2097];
    [neg_log_likeli_subj] = cal_log_likeli_optimal( best_params, data_subj, n_trials );
    log_likeli_subj = -neg_log_likeli_subj;
    
    log_likeli_opt = log_likeli_opt + log_likeli_subj;
    
    n_data_points_subj = size( data_subj, 1 );
    AIC_opt_subj(index_opt) = 2*(1) -2*log_likeli_subj;
    BIC_opt_subj(index_opt) = log(n_data_points_subj)*1 -2*log_likeli_subj;
    index_opt = index_opt + 1;
    
end

AIC_opt = 2*(subj_num) -2*log_likeli_opt;
BIC_opt = log(n_data_points)*subj_num -2*log_likeli_opt;

% cal likelihood of constant threshold model

log_likeli_const = 0;
params_subj0 = [100, 20];

for i = 1:subj_num
    
    data_subj = data_behavior( data_behavior(:,1) == i, : );
    best_params = fminsearch( @(params)...
        cal_log_likeli_const( params, data_subj, n_trials ), params_subj0 );
    
    % best_params = [120, 0.2077, 0.2097];
    [neg_log_likeli_subj] = cal_log_likeli_const( best_params, data_subj, n_trials );
    log_likeli_subj = -neg_log_likeli_subj;
    
    log_likeli_const = log_likeli_const + log_likeli_subj;
    
    n_data_points_subj = size( data_subj, 1 );
    AIC_const_subj(index_const) = 2*(2) -2*log_likeli_subj;
    BIC_const_subj(index_const) = log(n_data_points_subj)*2 -2*log_likeli_subj;
    index_const = index_const + 1;
    
end

AIC_const = 2*(subj_num*2) -2*log_likeli_const;
BIC_const = log(n_data_points)*subj_num*2 -2*log_likeli_const;

% cal likelihood of linear threshold model

log_likeli_linear = 0;
params_subj0 = [100, 0, 20];

for i = 1:subj_num
    
    data_subj = data_behavior( data_behavior(:,1) == i, : );
    best_params = fminsearch( @(params)...
        cal_log_likeli_linear( params, data_subj, n_trials ), params_subj0 );
    
    % best_params = [120, 0.2077, 0.2097];
    [neg_log_likeli_subj] = cal_log_likeli_linear( best_params, data_subj, n_trials );
    log_likeli_subj = -neg_log_likeli_subj;
    
    log_likeli_linear = log_likeli_linear + log_likeli_subj;
    
    n_data_points_subj = size( data_subj, 1 );
    AIC_linear_subj(index_linear) = 2*(3) -2*log_likeli_subj;
    BIC_linear_subj(index_linear) = log(n_data_points_subj)*3 -2*log_likeli_subj;
    index_linear = index_linear + 1;
    
end

AIC_linear = 2*(subj_num * 3) -2*log_likeli_linear;
BIC_linear = log(n_data_points)*subj_num*3 -2*log_likeli_linear;

% cal likelihood of optimal + linear threshold model

log_likeli_Olinear = 0;
params_subj0 = [0, 0, 20];

for i = 1:subj_num
    
    data_subj = data_behavior( data_behavior(:,1) == i, : );
    best_params = fminsearch( @(params)...
        cal_log_likeli_Olinear( params, data_subj, n_trials ), params_subj0 );
    
    % best_params = [120, 0.2077, 0.2097];
    [neg_log_likeli_subj] = cal_log_likeli_Olinear( best_params, data_subj, n_trials );
    log_likeli_subj = -neg_log_likeli_subj;
    
    log_likeli_Olinear = log_likeli_Olinear + log_likeli_subj;
    
    n_data_points_subj = size( data_subj, 1 );
    AIC_Olinear_subj(index_Olinear) = 2*(3) -2*log_likeli_subj;
    BIC_Olinear_subj(index_Olinear) = log(n_data_points_subj)*3 -2*log_likeli_subj;
    index_Olinear = index_Olinear + 1;
    
end

AIC_Olinear = 2*(subj_num * 3) -2*log_likeli_Olinear;
BIC_Olinear = log(n_data_points)*subj_num*3 -2*log_likeli_Olinear;

% cal likelihood of different thresholds model

log_likeli_thr = 0;
params_subj0 = 20;

for i = 1:subj_num
    
    data_subj = data_behavior( data_behavior(:,1) == i, : );
    thr_subj = thr( i, : );
    
    best_params = fminsearch( @(params)...
        cal_log_likeli_thr( params, data_subj, thr_subj, n_trials ), params_subj0 );
    
    [neg_log_likeli_subj] = cal_log_likeli_thr( best_params, data_subj, thr_subj, n_trials );
    log_likeli_subj = -neg_log_likeli_subj;
    
    log_likeli_thr = log_likeli_thr + log_likeli_subj;
    
    n_data_points_subj = size( data_subj, 1 );
    AIC_thr_subj(index_thr) = 2*(n_trials) -2*log_likeli_subj;
    BIC_thr_subj(index_thr) = log(n_data_points_subj)*n_trials -2*log_likeli_subj;
    index_thr = index_thr + 1;
    
end

AIC_thr = 2*(subj_num*n_trials) -2*log_likeli_thr;
BIC_thr = log(n_data_points)*subj_num*n_trials -2*log_likeli_thr;
% cal likelihood of optimal+C model

log_likeli_opt_C = 0;
params_subj0 = [0, 20];

for i = 1:subj_num
    
    data_subj = data_behavior( data_behavior(:,1) == i, : );
    
    best_params = fminsearch( @(params)...
        cal_log_likeli_optimal_C( params, data_subj, n_trials ), params_subj0 );
    
    % best_params = [120, 0.2077, 0.2097];
    [neg_log_likeli_subj] = cal_log_likeli_optimal_C( best_params, data_subj, n_trials );
    log_likeli_subj = -neg_log_likeli_subj;
    
    log_likeli_opt_C = log_likeli_opt_C + log_likeli_subj;
    best_params_all_OC(i, :) = best_params;
    
    n_data_points_subj = size( data_subj, 1 );
    AIC_OC_subj(index_OC) = 2*(2) -2*log_likeli_subj;
    BIC_OC_subj(index_OC) = log(n_data_points_subj)*2 -2*log_likeli_subj;
    index_OC = index_OC + 1;
    
end

AIC_OC = 2*( subj_num*2 ) - 2*log_likeli_opt_C;
BIC_OC = log(n_data_points)*subj_num*2 -2*log_likeli_opt_C;

end
subj_num

% LR_CT_opt = -2*( log_likeli_opt - log_likeli_const )
% p1 = 1 - cdf( 'chi2', LR_CT_opt, subj_num )
% 
% LR_3p_opt = -2*( log_likeli_opt - log_likeli_3p )
% p2 = 1 - cdf( 'chi2', LR_3p_opt, 2*subj_num )
% 
% LR_4p_opt = -2*( log_likeli_opt - log_likeli_4p )
% p3 = 1 - cdf( 'chi2', LR_4p_opt, 3*subj_num )
% 
% LR_3p_const = -2*( log_likeli_const - log_likeli_3p )
% p4 = 1 - cdf( 'chi2', LR_3p_const, subj_num )
% 
% LR_4p_const = -2*( log_likeli_const - log_likeli_4p )
% p5 = 1 - cdf( 'chi2', LR_4p_const, 2*subj_num )
% 
% LR_4p_3p = -2*( log_likeli_3p - log_likeli_4p )
% p6 = 1 - cdf( 'chi2', LR_4p_3p, subj_num )
% 
% LR_thr_4p = -2*( log_likeli_4p - log_likeli_thr )
% p7 = 1 - cdf( 'chi2', LR_thr_4p, (n_trials-4)*subj_num )


t_dist_multi = icdf( 't', 0.975, subj_num - 1 );

AIC_3p_4p_subj = AIC_3p_subj - AIC_4p_subj;
AIC_opt_4p_subj = AIC_opt_subj - AIC_4p_subj;
AIC_OC_4p_subj = AIC_OC_subj - AIC_4p_subj;
AIC_const_4p_subj = AIC_const_subj - AIC_4p_subj;
AIC_thr_4p_subj = AIC_thr_subj - AIC_4p_subj;
AIC_linear_4p_subj = AIC_linear_subj - AIC_4p_subj;
AIC_Olinear_4p_subj = AIC_Olinear_subj - AIC_4p_subj;


% AIC_all = [ AIC_3p_subj', AIC_4p_subj', AIC_opt_subj', AIC_const_subj',...
%     AIC_OC_subj', AIC_thr_subj', AIC_linear_subj', AIC_Olinear_subj' ];

% AIC_delta_all = [ AIC_3p_4p_subj', AIC_thr_4p_subj', AIC_OC_4p_subj', AIC_const_4p_subj',...
%     AIC_opt_4p_subj', AIC_linear_4p_subj', AIC_Olinear_4p_subj' ]; 

AIC_delta_all = [ AIC_3p_4p_subj', AIC_thr_4p_subj', AIC_linear_4p_subj', AIC_const_4p_subj'...
    , AIC_OC_4p_subj', AIC_opt_4p_subj']; 

AIC_all = [ AIC_4p_subj', AIC_3p_subj', AIC_thr_subj', AIC_linear_subj'...
    , AIC_const_subj',AIC_OC_subj', AIC_opt_subj'];



figure1 = figure('PaperUnits','centimeters','PaperSize',[13 13],...
    'PaperPosition',[0 0 13 13],'Units','centimeters',...
    'Position',[0 0 13 13],'Color',[1 1 1]);

axes1 = axes('Parent',figure1);
hold(axes1,'on');
ylim([-50, 400])

violinplot( AIC_delta_all )
plot( [0.5, 6.5], [0, 0],'k' )


for i = 1:6
   
    [CI] = cal_CI_boot( AIC_delta_all(:,i), 100000 );
    
    AIC_delta_95CI(1,i) = CI(1);
    AIC_delta_95CI(2,i) = CI(3);
    
end

AIC_delta_95CI

    





