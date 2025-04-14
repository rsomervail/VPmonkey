%
%  script to split raw bdf files with dual EEG into two EEGlab files that can be loaded by the main script 
%  
%%
clc
clearvars
% addlab
addpath '/home/rick/Dropbox/Somervail_Richard/VPmonkey/scripts';

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

% homedir = '/home/rick/neuro/iannettilab/Projects/VP monkey/data/raw_dual2split';
homedir = '/media/rick/Rick_LabData_3/neuro/iannettilab_ongoing/VPmonkey/data/raw_dual2split';
cd(homedir);

s.savePath =  ['/media/rick/Rick_LabData_3/neuro/iannettilab_ongoing/VPmonkey/data/raw_' subfold]; % mkdir(s.savePath)

s.loadEyetrack = true;
s.loadLFP      = true; if s.loadLFP, addpath '/home/rick/Dropbox/Somervail_Richard/VPmonkey/scripts/TDTMatlabSDK/TDTSDK/TDTbin2mat'; end
s.loadSPK      = true;

s.conds = {'VIS', 'SOM', 'AUD'};
s.lfp_trig      = [2, 8, 32];
s.lfp_freqs     = [1, 100]; 
s.spike_freqs   = [100, 2000]; 

s.chans2remove  = { 'EXG4','EXG5','EXG6','EXG7','EXG8', 'X1','X2','X3','X4'  }; % pilot data
% s.chans2remove  = { 'EXG4','EXG5','EXG6','EXG7','EXG8', 'X1','X2','X3','X4'  }; % main experiment

%% select sessions
sesh = dir; sesh = {sesh.name}; sesh(1:2) = [];
sel = listdlg('ListString', sesh, 'SelectionMode','multiple');
sesh = sesh(sel);
nsesh = length(sesh);

%% loop through sessions
tIN_sesh = tic;
for sh = 1 :nsesh

    % grab raw data file names
    cd([homedir filesep sesh{sh}])
    files = dir; files = {files.name}; files(~endsWith(files, '.bdf')) = [];
    nfiles = length(files);

    numChans = 94; % initial number of channels to load, can change after loading one file

    % grab matching eye tracking files
    if s.loadEyetrack 
        files_eye = dir; files_eye(~startsWith({files_eye.name}, '2022')) = [];
        [~, neworder] = sort([files_eye.datenum]);
        files_eye = files_eye(neworder); clear neworder
        files_eye = {files_eye.name};
    end

    %% loop through files (blocks)
    SPKSORTm = cell(nfiles,1);
    SPKSORTt = cell(nfiles,1);
    for f = 1 :nfiles

        % load EEG file 
        try 
            EEG = pop_biosig( files{f} , 'importevent','off', 'rmeventchan' ,'off', 'channels', 1:numChans);
        catch
            disp 'checking number of channels'
            EEG = pop_biosig( files{f} , 'importevent','off'); % check how many channels there are 
            numChans = EEG.nbchan + 1; % get number of channels (excluding status)
            EEG = pop_biosig( files{f} , 'importevent','off', 'rmeventchan' ,'off', 'channels', 1:numChans); 
        end

        % extract triggers
        if any(strcmp({EEG.chanlocs.labels}, 'Status'))
            trigChan = EEG.data(strcmp({EEG.chanlocs.labels}, 'Status'),:);
        else
            error 'status (trigger) channel not loaded - perhaps load function above is not loading all channels'
        end
        trigChanDiff = diff(trigChan); %
        triggers  = find(trigChanDiff > 0) + 1; % find timepoints where the trigger channel increases
        trigCodes = trigChanDiff( triggers -1 ); % find the amount the trigger channel changes by (this is the code of the trigger)
        %     figure; plot(trigChan-mode(trigChan)); hold on; ylim([-1 6]); xlim([1802262-1000, 1802262+1000])
        trigs2remove = trigCodes > 10; % RS: some triggers are clearly artifacts because they are so big, so I'm removing any triggers which are above 10
        triggers(trigs2remove) = [];
        trigCodes(trigs2remove) = [];
        clear trigs2remove
        trigCodes(  trigCodes > 3  ) = 3; % RS: one of the auditory triggers was the wrong value, so I'm correcting them all to "3" here
        while true
            if length(trigCodes) > 30
                trigCodes(1) = [];
                triggers(1)  = [];
            else 
                break
            end
        end
        % loop through triggers and construct events structure
        events = struct;
        for trl = 1 : length(triggers)
            events(trl).type = s.conds{trigCodes(trl)};
            events(trl).urevent =  trl;
            events(trl).latency =  triggers(trl);
        end
        EEG.event = events; clear events % add events to EEG structure
        EEG.urevent = EEG.event; EEG.urevent = rmfield(EEG.urevent, 'urevent');
        EEG = eeg_checkset(EEG); 

        % remove unused EEG channels (but leave missing channels we want to reconstruct)
        EEG = pop_select( EEG, 'nochannel',  find( cellfun( @(x) contains(x, s.chans2remove) , {EEG.chanlocs.labels} ) ) );

        % split data into two
        chansM = find(startsWith({EEG.chanlocs.labels},'1-') & strcmp({EEG.chanlocs.type},'EEG')); % get Miguele EEG chans
        chansT = find(startsWith({EEG.chanlocs.labels},'2-') & strcmp({EEG.chanlocs.type},'EEG')); % get Tulio EEG chans
        chanPhoto = find(strcmp({EEG.chanlocs.labels},'1-Erg1')); EEG.chanlocs(chanPhoto).labels = 'photo'; EEG.chanlocs(chanPhoto).type = 'Misc';
        chanMicro = find(strcmp({EEG.chanlocs.labels},'2-Erg1')); EEG.chanlocs(chanMicro).labels = 'micro'; EEG.chanlocs(chanMicro).type = 'Misc';
        EEGm = pop_select(EEG, 'channel', [chansM chanPhoto chanMicro] ); % Miguele
        EEGt = pop_select(EEG, 'channel', [chansT chanPhoto chanMicro] ); % Tulio
        clear EEG

        % remove 1- and 2- from the chanlabels      ctemp.
        for c = 1:EEGm.nbchan, EEGm.chanlocs(c).labels = strrep( EEGm.chanlocs(c).labels, '1-','' ); end
        for c = 1:EEGt.nbchan, EEGt.chanlocs(c).labels = strrep( EEGt.chanlocs(c).labels, '2-','' ); end

        % reorder so that EEG electrodes come first
        newOrderT =  [ find(strcmp({EEGt.chanlocs.type},'EEG')) , find(strcmp({EEGt.chanlocs.type},'Misc')) ];
        EEGt.chanlocs = EEGt.chanlocs(newOrderT);
        EEGt.data = EEGt.data(newOrderT,:);

        % load eye tracking data and add channels to EEG data
        if s.loadEyetrack
            disp 'loading & aligning eye tracker data ...'
            [eyeData, eyeEvents] = VPmonkey_importEyetrackingData(files_eye{f});
            eyeEvents = eyeEvents(ismember(  eyeEvents.DeltaTime,  s.conds   ) , :);
            
            % split by subject
            indym = endsWith( eyeData.Properties.VariableNames , '1');
            indym(1) = false; indym(end-1) = true;
            indyt = ~indym; indyt([1,end]) = false; indyt(end-1) = true;
            eyeDatam = eyeData(:,indym);
            eyeDatat = eyeData(:,indyt);
            disp '... loaded eyetracker'
          
            % merge data with EEG (aligning, interpolating, cropping and concatenating)
            disp 'appending eyetracker channels to EEG data ...'
            EEGm = VPmonkey_combineEEGeye(EEGm, eyeDatam, eyeEvents);
            EEGt = VPmonkey_combineEEGeye(EEGt, eyeDatat, eyeEvents);
            disp '... appended eyetracker'

        end

        % load intracranial data and extract LFPs, then add channels to EEG data
        if s.loadLFP
            disp 'loading LFP ...'
            lfp_path = [cd filesep files{f}(1:end-4) ];
            [LFPm, LFPt] = VPmonkey_loadLFP( lfp_path  , s );
            clear lfp_path
            disp '... loaded LFP'

            % rename channels
            for c = 1:5
                LFPm.chanlocs(c).labels = ['LFP_' 'el' num2str(c) ];
                LFPt.chanlocs(c).labels = ['LFP_' 'el' num2str(c) ];
%                 LFPm.chanlocs(c).labels = ['LFP_' sesh{sh}(1:3) '_el' num2str(c) ];
%                 LFPt.chanlocs(c).labels = ['LFP_' sesh{sh}(1:3) '_el' num2str(c) ];
                LFPm.chanlocs(c).type   = 'LFP';
                LFPt.chanlocs(c).type   = 'LFP';
            end

            disp 'aligning and appending LFP channels to EEG data ...'
            EEGm = VPmonkey_combineEEGlfp(EEGm, LFPm);
            EEGt = VPmonkey_combineEEGlfp(EEGt, LFPt);
            disp '... appended LFP'
        end

        % save EEG/LFP/EYE
        savePathM = [s.savePath filesep  sesh{sh} ' SubM' filesep]; mkdir(savePathM)
        savePathT = [s.savePath filesep  sesh{sh} ' SubT' filesep]; mkdir(savePathT)
        saveNameM = [ files{f}(1:end-4) ' SubM' ]; EEGm.setname = saveNameM;
        saveNameT = [ files{f}(1:end-4) ' SubT' ]; EEGt.setname = saveNameT;
        pop_saveset(EEGm, 'filepath', savePathM, 'filename', [saveNameM ] );
        fprintf('... saved session %s, file %s (Miguele): %s\n', sesh{sh}, files{f}, [savePathM saveNameM]) 
        pop_saveset(EEGt, 'filepath', savePathT, 'filename', [saveNameT ]);
        fprintf('... saved session %s, file %s (Tulio): %s\n',   sesh{sh}, files{f}, [savePathT saveNameT])

        % LOAD SPIKES  
        if s.loadSPK
            disp 'loading spikes ...'
            lfp_path = [cd filesep files{f}(1:end-4) ];
            [SPKm, SPKt] = VPmonkey_loadSPK( lfp_path  , s );
            clear lfp_path
            disp '... loaded spikes'

            % rename channels
            for c = 1:5
                SPKm.chanlocs(c).labels = ['SPK_' 'el' num2str(c) ];
                SPKt.chanlocs(c).labels = ['SPK_' 'el' num2str(c) ];
                SPKm.chanlocs(c).type   = 'SPK';
                SPKt.chanlocs(c).type   = 'SPK';
            end

            % resample spikes to round number 
            SPKm = pop_resample(SPKm, 24000);
%             SPKm = pop_rs_downsample(SPKm, 2); % no longer downsampling to make spike detection more accurate
            SPKt = pop_resample(SPKt, 24000);
%             SPKt = pop_rs_downsample(SPKt, 2); 
            % round event latencies after downsampling
            for e = 1:length(SPKm.event)
                SPKm.event(e).latency = round(SPKm.event(e).latency);
            end
            for e = 1:length(SPKt.event)
                SPKt.event(e).latency = round(SPKt.event(e).latency);
            end

            % align with dummy EEG (for the events)
            disp 'aligning SPK channels to EEG ...'
            tempm = pop_select(EEGm,'channel',find(strcmp({EEGm.chanlocs.labels},'CZ')));
            tempm = pop_resample(tempm, SPKm.srate);
            tempt = pop_select(EEGt,'channel',find(strcmp({EEGt.chanlocs.labels},'CZ')));
            tempt = pop_resample(tempt, SPKt.srate);
            SPKm = VPmonkey_combineEEGlfp(tempm, SPKm);
            SPKt = VPmonkey_combineEEGlfp(tempt, SPKt);
            SPKm = pop_select(SPKm,'channel',find(~strcmp({SPKm.chanlocs.labels},'CZ')));
            SPKt = pop_select(SPKt,'channel',find(~strcmp({SPKt.chanlocs.labels},'CZ')));

            disp '... aligned spikes'

            % SAVE
            savePathM = [s.savePath filesep  sesh{sh} ' SubM' filesep]; mkdir(savePathM)
            savePathT = [s.savePath filesep  sesh{sh} ' SubT' filesep]; mkdir(savePathT)
            saveNameM = [ 'SPK ' files{f}(1:end-4) ' SubM' ]; SPKm.setname = saveNameM;
            saveNameT = [ 'SPK ' files{f}(1:end-4) ' SubT' ]; SPKt.setname = saveNameT;
            pop_saveset(SPKm, 'filepath', savePathM, 'filename', [ saveNameM ] );
            fprintf('... saved session %s, file %s (Miguele): %s\n', sesh{sh}, files{f}, [savePathM saveNameM]) 
            pop_saveset(SPKt, 'filepath', savePathT, 'filename', [ saveNameT ]);
            fprintf('... saved session %s, file %s (Tulio): %s\n',   sesh{sh}, files{f}, [savePathT saveNameT])

            % store SPK data for concatenation and export for spike sorting
            SPKSORTm{f} = SPKm;
            SPKSORTt{f} = SPKt;

        end

        fprintf('finished file %d/%d (tElapsed = %.2f mins)\n', f, nfiles, toc(tIN_sesh)/60  )
    end % end file (block) loop

    %% concatenate blocks for export to spike sorting
    if s.loadSPK
        
        % concatenate
        mdata   = [];
        mtimes  = [];
        tdata   = [];
        ttimes  = [];
        for f = 1:nfiles
            % sub M
            mdata   = [mdata, SPKSORTm{f}.data  ];
            mtimes  = [mtimes, SPKSORTm{f}.times/1000 ];
            % sub T
            tdata   = [tdata, SPKSORTt{f}.data  ];
            ttimes  = [ttimes, SPKSORTt{f}.times/1000 ];

        end % file loop

        % export
        save([savePathM filesep 'SPKSORT data ' sesh{sh} ' SubM'], 'mdata')
        save([savePathT filesep 'SPKSORT data ' sesh{sh} ' SubT'], 'tdata')
        save([savePathM filesep 'SPKSORT times ' sesh{sh} ' SubM'], 'mtimes')
        save([savePathT filesep 'SPKSORT times ' sesh{sh} ' SubT'], 'ttimes')

    end
end % end session loop

cd([ s.savePath   ])
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
