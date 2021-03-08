function data_behavior = format_data( data_all )

subj_num = size( data_all, 1 );
round_num = size( data_all, 2 );

n_trials = 0;
for subjNum = 1:subj_num
    for rnd = 1:round_num
        n_trials = n_trials + data_all{subjNum, rnd}.trialNumber;
        % n_trials = n_trials + 1;
    end
end


running_index = 1;
data_behavior = NaN( n_trials, 5 );

for subjNum = 1:subj_num
    for rnd = 1:round_num
        
        for trialIndex = 1:data_all{subjNum, rnd}.trialNumber
            % for trialIndex = 1
            value =  data_all{subjNum, rnd}.value(trialIndex);
            took_it = trialIndex == data_all{subjNum, rnd}.trialNumber;
            data_behavior( running_index, : ) = [ subjNum, rnd, trialIndex, value, took_it ];
            
            running_index = running_index + 1;
            
        end
    end
end
