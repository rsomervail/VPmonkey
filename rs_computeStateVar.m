


% function out = rs_computeStateVar( states , win )

% states = EEG.microstate.labels;
win = round(   ( 400/1000) * 2048 );

    difs = diff( states, [], 2 ) ~= 0;

    out = nan(size(states));
%     for trl = 1:size(states,1)
% 

        for t = 1:size(states,2)

            try 
                out(:,t) = sum( difs(:,  t-win : t+win) , 2 );
                t
            end
        
        end

%     end




% end

% figure; plot(EEG.times,out)