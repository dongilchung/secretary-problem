function [CI] = cal_CI_boot( samples, n_samples )


if size( samples, 1 ) < size( samples, 2 )
    samples = samples';
end
    
sample_size = length( samples );
sm = zeros( n_samples, 1 );

for i = 1:n_samples
    
    ri = randi( sample_size, 1, sample_size );
    si = samples( ri );
    
    sm(i) = mean( si );
    
end

sorted_sm = sort( sm );

CI_low = sorted_sm( round( 0.025 * n_samples ) );
CI_high = sorted_sm( round( 0.975 * n_samples ) );

CI = [CI_low, mean(samples), CI_high];

