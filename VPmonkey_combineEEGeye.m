% 
%
%    
% % EEGlab data, eye tracker data, eye tracker events
%
%%
function EEG = VPmonkey_combineEEGeye(EEG, eye, ents)

% get sampling rate
fs = EEG.srate; 

% rename variables which end in a 1
if ismember('TotalTime1', eye.Properties.VariableNames)
    eye = renamevars(eye, eye.Properties.VariableNames  , strrep(eye.Properties.VariableNames,'1',''));
end

% get time vector for eye tracker 
start_indy = find([eye.Count]==1);
eye = eye(start_indy:end,:);
times_eye  = eye.TotalTime;

% remove nans/infs
indy2remove = ~isfinite(times_eye);
times_eye(indy2remove) = [];
eye(indy2remove,:) = [];

% check for repeat timestamps and remove them
repeatVals = str2double(eye.DeltaTime)==0;
if any(repeatVals)
    eye(repeatVals,:) = [];
    times_eye(repeatVals) = [];
end

% check that EEG timecourse is longer than eye timecourse (it should be from the protocol we used)
times_eeg = EEG.times/1000;
% if times_eeg(end) < times_eye(end), error 'RS error: time mismatch between eye and eeg'; end 
if times_eeg(end) < times_eye(end)
    truncateEye = true; 
else
    truncateEye = false; 
end 

% % convert time vectors to fractions of milliseconds
% df = 10 * 1000;
% times_eye  = round(   times_eye *df   );
% times_eeg  = round(   times_eeg'*df    ); 

%% standard interpolation functions 

% upsample eye time vector
times_eye_up = [times_eye(1) : 1/fs : times_eye(end)]'; 
% times_eye_up = [times_eye(1) : 1/(fs/df) : times_eye(end)]';  % for use with df above, changing units
% times_eye_up = linspace( times_eye(1) , times_eye(end), ceil( (times_eye(end) - times_eye(1)) * fs ) )' ; % upsample eye time vector

% linearly interpolate each signal to match srate & make regular timing
flds = eye.Properties.VariableNames; flds( contains(flds,'Time') | strcmp(flds,'Count') ) = []; 
eye_up = table(times_eye_up, 'VariableNames', {'time'}); 
for f = 1  :length(flds)
    temp = eye.(flds{f});

    temp_up = interp1( times_eye, temp, times_eye_up, 'nearest'); 
    
    eye_up.(flds{f}) = temp_up;
end


%% figure out how many values to pad with
% check that events match
ents_eeg = EEG.event;
if ~isequaln( {ents_eeg.type}', ents.DeltaTime  )
    error 'events do not match!!!!'
else
    nents = length(EEG.event);
end

% get EEGlab event times in seconds
ents = table2struct(ents);
for e = 1:nents
    ents_eeg(e).time = EEG.times(ents_eeg(e).latency)/1000;
end

% compute mean difference in times between corresponding events & use to calculate number of preceding samples necessary
diffs = [ents_eeg.time] - [ents.TotalTime];
leadingTime = mean(diffs);
numLeadingSamps = round(mean(diffs) * fs);

% if necessary, truncate eye tracking data (in the case of prematurely ended blocks)
if truncateEye

    overshootSamps = (length(times_eye_up) + numLeadingSamps) - EEG.pnts;
    eye_up( end-overshootSamps+1 : end, : ) = [];

%     ? old method:
%     [~, lastIndy] = min(abs( times_eeg(end) - times_eye_up ));
%     eye_up = eye_up(1:lastIndy,:);
% %     times_eye_up = times_eye_up(1:lastIndy); % not necessary to do this unless order of operations changes

    numTrailingSamps = 0;
else % normal block
    % compute lagging samples (if this was not a block where the eye data outlasted the eeg data)
    numTrailingSamps = EEG.pnts - (length(times_eye_up) + numLeadingSamps);
end

%% combine pupil height and width into pupil diameter
eye_up.pupilDiameter = mean(  [eye_up.PupilHeight eye_up.PupilWidth ], 2 );
eye_up = eye_up( :, [1:8 11 9:10] ); % reorder

%% add eye channels to EEG data
flds = eye_up.Properties.VariableNames; flds = flds(2:end);
for f = 1:length(flds)

    temp = eye_up.(flds{f});
    temp = [ nan(1,numLeadingSamps), temp', nan(1,numTrailingSamps)   ]; % pad missing data with nans
    
    % add to EEG
    EEG.data(end+1,:) = temp;
    EEG.chanlocs(end+1).labels = ['eye_' flds{f}]; 
    EEG.chanlocs(end).type   = 'Eye'; 

end
EEG.nbchan = EEG.nbchan + length(flds); % update number of channels



end