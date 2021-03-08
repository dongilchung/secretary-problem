clear
clc

% load('data_all.mat')
load('data_all_exp4.mat')
subj_num = size( data_all, 1 );
% data_behavior = format_data( data_all );
clear data_all
% load('thr_exp3_22.mat')
% thr = thr_22;

load('UBSDM_Exp4_p3_r1.mat')
n_params = 3;
% load('thr.mat')

grid_length = 40;

for subjIndex = 1:subj_num
    % for paramIndex = 1:n_params
        subjIndex
      %  runningIndex = (subjIndex-1)*n_params + paramIndex;
        params1Index = subjIndex*n_params - 2;
        params3Index = subjIndex*n_params;
        doi = samples( params1Index:params3Index, 2001:end )';
        xGrid = zeros( grid_length, n_params );
        
        for i = 1:n_params
            xGrid(:,i) = linspace( min(doi(:,i)), max(doi(:,i)), grid_length );
        end
        
        [x1, x2, x3] = ndgrid( xGrid(:,1), xGrid(:,2), xGrid(:,3) ); 
        xi = [x1(:), x2(:), x3(:)];
        sigma4 = std( doi );
        
        bw = sigma4*(4/((n_params+2)*length(doi)))^(1/(n_params+4));
        
        f = mvksdensity( doi, xi, 'bandwidth', bw );
        
        % [d,i] = ksdensity( doi, 'NumPoints', 300 );
        
        [~,maxI] = max(f);
        best_params_all( subjIndex, : ) = xi(maxI,:);
    % end
    
%     data_subj = data_behavior( data_behavior(:,1) == subjIndex, : );
%     
%     [log_likeli_subj, threshold, exp_value] = cal_log_likeli_subj_wc( best_params_all(subjIndex,:), data_subj );
%     
%     thresh_all(subjIndex,:) = threshold;
%     exp_value_all(subjIndex) = exp_value;
    
end

% thr_0 = [thresh_all, zeros(subj_num,1)];
%
%
% close all
% for i = 1:subj_num
% 
%     subplot( 4,6,i ),plot( thr_0(i,:), '-or' )
%     hold on
%     subplot( 4,6,i ),plot( thr(i,:), '-ob' )
%     subplot( 4,6,i ),plot( [0,5],[best_params_all(i,1), best_params_all(i,1)], '-k' )
%     ylim([0,150])
% 
% end

best_40_Exp4_p3 = best_params_all;
save best_40_Exp4_p3_r1 best_40_Exp4_p3



%
% load('big5_result.mat')
% load('PAL_LAB_lotr.mat')
% [r,p] = corr( thr_all, five_scores )
% [r,p] = corr( best_params_all, five_scores )
% [r,p] = corr( best_params_all, lotr_total )
% [r,p] = corr( thr_all, lotr_total )
%
%
% load('earning.mat')
% [r,p] = corr( thr_all(:,1)-thr_all(:,4), earn' )
%
%






