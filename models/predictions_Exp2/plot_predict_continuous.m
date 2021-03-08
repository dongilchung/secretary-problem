clear
clc
close all
load('prediction_results_p4.mat')
load('thresh_exp2.mat')
load('thresh_model_exp2_5')
load('continuous_rp.mat')


[~, exp_value5] = toy_model_power( 1, 100, 5 );
[~, exp_value2] = toy_model_power( 1, 100, 2 );
[~, exp_value10] = toy_model_power( 1, 100, 10 );

x_2_5 = exp_value2 - exp_value5;
x_10_5 = exp_value10 - exp_value5;



max_boundary_temp = max_boundary;
min_boundary_temp = min_boundary;

for i = 1:length(max_boundary)
    
    index1 = max( 1, i-1);
    index2 = min( length(max_boundary), i + 1 );
    max_boundary(i) = mean( max_boundary_temp( index1:index2 ) );
    min_boundary(i) = mean( min_boundary_temp( index1:index2 ) );
    
end


exp_diff = exp_diff(2:end-1);
max_boundary = max_boundary(2:end-1);
min_boundary = min_boundary(2:end-1);


figure1 = figure('PaperUnits','centimeters','PaperSize',[10 13],...
    'PaperPosition',[0 0 10 13],'Units','centimeters',...
    'Position',[0 0 10 13],'Color',[1 1 1]);

axes1 = axes('Parent',figure1);
hold(axes1,'on');
ylim([65, 105])

mean_predict = 0.5*( max_boundary + min_boundary );
confplot( exp_diff, mean_predict, max_boundary-mean_predict, mean_predict - min_boundary  );


% 
% 
predict_mean_o2 = (max_boundary + min_boundary)/2;
predict_ci_o2 = max_boundary - (max_boundary + min_boundary)/2;

n_subj_2 = size( thr2, 1 );
n_subj_5 = size( thr5, 1 );
last_1_2 = thr2(:,1);
last_1_5 = thr5(:,4);


load('prediction10.mat')
load('thr_exp4.mat')
thr10 = thr;

predict_mean_o10 = (max_boundary + min_boundary)/2;
predict_ci_o10 = max_boundary - (max_boundary + min_boundary)/2;

diff_9th = - mean(mean( thr10( thr10(:,9) > 0, 7:8 ), 2 ) - thr10( thr10(:,9) > 0, 9 ) )
thr10( thr10(:,9) < 0, 9) = mean( thr10( thr10(:,9) < 0, 7:8 ), 2 ) + diff_9th;

n_subj_10 = size( thr10, 1 );
last_1_10 = thr10(:,9);

% 
% 
% figure1 = figure('PaperUnits','centimeters','PaperSize',[6 13],...
%     'PaperPosition',[0 0 6 13],'Units','centimeters',...
%     'Position',[0 0 6 13],'Color',[1 1 1]);
% 
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');


% 
% for i = 1:n_subj_2
%     plot( 1, last_1_2(i), 'k-o','MarkerSize',6,'linewidth',1 );
%     hold on
% end
% 
% for i = 1:n_subj_5
%     plot( 2, last_1_5(i), 'k-o','MarkerSize',6,'linewidth',1 );
%     hold on
% end
% 
% for i = 1:n_subj_10
%     plot( 3, last_1_10(i), 'k-o','MarkerSize',6,'linewidth',1 );
%     hold on
% end

% 
% 
% ylim([50, 150])
% xlim([0.5, 3.5])

% Y = [ thr2(:,1); thr5(:,4); thr10(:,9)];
% G = [ ones(length(thr2),1); 2*ones(length(thr5),1); 3*ones(length(thr10),1)];
% 
% [p,tbl,stats] = anova1( Y, G, 'off' )
% COMPARISON = multcompare( stats )

hold on

Y2 = icdf('t', 0.025, n_subj_2-1);
Y5 = icdf('t', 0.025, n_subj_5-1);
Y10 = icdf('t', 0.025, n_subj_10-1);

mean_obser = [ mean( last_1_2 ), mean( last_1_5 ), mean( last_1_10 ) ];
errorbar_obser = [ Y2*std(last_1_2)/sqrt(n_subj_2), Y5*std(last_1_5)/sqrt(n_subj_5),...
    Y10*std(last_1_10)/sqrt(n_subj_10)];


h = errorbar( [x_2_5, 0, x_10_5], mean_obser, errorbar_obser,...
    'ro', 'markerfacecolor','w', 'linewidth',2,'MarkerSize',8);
h.CapSize = 0;

% 
% 
% hold on
% 
% estimate_o5_mean = mean( thresh_all(:,4) );
% estimate_o5_ci = Y5*std( thresh_all(:,4) )/sqrt(n_subj_5);
% 
% mean_predict = [ predict_mean_o2, estimate_o5_mean, predict_mean_o10 ];
% ci_predict = [ predict_ci_o2, estimate_o5_ci, predict_ci_o10 ];
% 
% 
% % mean_thr = mean( last_1 );
% 
% h = errorbar( [1.9], mean_predict(2), ci_predict(2), 'go',...
%     'markerfacecolor','w', 'linewidth',2,'MarkerSize',8);
% h.CapSize = 0;
% 
% h = errorbar( [0.9, 2.9], mean_predict([1,3]), ci_predict([1,3]), 'bo',...
%     'markerfacecolor','w', 'linewidth',2,'MarkerSize',8);
% h.CapSize = 0;







