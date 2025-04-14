%
%
%
%
%
%
function EEG = VPmonkey_combineEEGlfp(EEG, LFP)

% check that events are the same
if ~isequaln({EEG.event.type},{LFP.event.type})
    error 'RS ERROR: EVENT MISMATCH BETWEEN EEG AND LFP, NEED EXCEPTIONS'
end

% compute median time-difference (in sample units) between matching events
lats_eeg = [EEG.event.latency];
lats_lfp = [LFP.event.latency];
timediff = lats_eeg - lats_lfp;
mdiff = round(median(timediff));

% insert LFP chans to EEG to align  
temp2add = [];
if      mdiff > 0  % positive means eeg events are later than their LFP counterparts, so pad LFP with mdiff nans 
    temp2add = nan(5,mdiff); % pad start
    endPointDiff = (LFP.pnts + mdiff)-EEG.pnts;
    if      endPointDiff > 0 % if after padding the LFP at the start, the LFP is now too long for the EEG, then truncate at end
        temp2add = [ temp2add  LFP.data(:,1:end-endPointDiff)  ];
    elseif  endPointDiff < 0 % if after padding the LFP, the LFP is still too short for the EEG then pad at the end too
        temp2add = [ temp2add  LFP.data  nan(5, abs(endPointDiff))  ];
    else % otherwise just add the whole data in with no padding at the end
        temp2add = [ temp2add  LFP.data ];
    end

elseif  mdiff < 0 % negative means eeg events are earlier than the LFP events, so truncate start of LFP
    
    endPointDiff = (LFP.pnts + mdiff)-EEG.pnts;
    if      endPointDiff > 0 % if after truncating the LFP at the start, the LFP is still too long for the EEG, then truncate at end
        temp2add = LFP.data(:, abs(mdiff)+1 : end-endPointDiff); 
    elseif  endPointDiff < 0 % if after truncating the LFP, the LFP is now too short for the EEG then pad at the end too
        temp2add = [ LFP.data(:, abs(mdiff)+1 :end)   nan(5, abs(endPointDiff))  ]; 
    else % otherwise just add the whole data in with no padding at the end
        temp2add =   LFP.data(:, abs(mdiff)+1 :end);  
    end


end
EEG.data = [EEG.data; temp2add];

% set channel labels and type
for c = 1:LFP.nbchan
    EEG.chanlocs(EEG.nbchan+c).labels = LFP.chanlocs(c).labels;
    EEG.chanlocs(EEG.nbchan+c).type   = LFP.chanlocs(c).type;
end
EEG.nbchan = EEG.nbchan + LFP.nbchan;

end % end function
