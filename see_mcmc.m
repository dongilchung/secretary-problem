clear
clc


load('data_all_exp2_10.mat')
subj_num = size( data_all, 1 );
data_behavior = format_data( data_all );
clear data_all

load('UBSDM_Exp4_r1.mat')
n_params = 4;

grid_length = 40;

for subjIndex = 1:subj_num
    % for paramIndex = 1:n_params
        subjIndex
      %  runningIndex = (subjIndex-1)*n_params + paramIndex;
        params1Index = subjIndex*n_params - 3;
        params4Index = subjIndex*n_params;
        doi = samples( params1Index:params4Index, 2001:end )';
        xGrid = zeros( grid_length, n_params );
        
        for i = 1:n_params
            xGrid(:,i) = linspace( min(doi(:,i)), max(doi(:,i)), grid_length );
        end
        
        [x1, x2, x3, x4] = ndgrid( xGrid(:,1), xGrid(:,2), xGrid(:,3), xGrid(:,4)); 
        xi = [x1(:), x2(:), x3(:), x4(:)];
        sigma4 = std( doi );
        bw = sigma4*(4/(length(doi)+2))^(1/(length(doi)+4));
        
        f = mvksdensity( doi, xi, 'bandwidth', bw );
        
        % [d,i] = ksdensity( doi, 'NumPoints', 300 );
        
        [~,maxI] = max(f);
        best_params_all( subjIndex, : ) = xi(maxI,:);
    
end


best_40_Exp4_3 = best_params_all;
save best_40_Exp4 best_40_Exp4


