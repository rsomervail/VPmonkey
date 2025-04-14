%
%  
%  ! script to load and split the two bdf files by channels somehow, then rename channels when loading here
%  ! correct onset delay of each stimulus modality using the two Erg channels and electrical stimulus artifact
%       ! don't correct trial-by-trial jitter as this seems small but definitely estimate it and save the trial-by-trial values & estimate for possible later use
%  ! make sure that channels are interpolated similarly to how this is done in Giacomo's thing (maybe even do within fieldtrip using his functions)
%       ! right now Cz is being interpolated weird, maybe including some far away edge electrodes, seems like a gap relative to Fcz which is huge ERP
%  ! manually check interpolated channels (plot a subplot for each chan with auto limits, plot all 28 electrodes mark bad ones in red plots, do 30 = 5*6 subplots )
%       ! edge electrodes on the top right seem to still be bad after interpolation - why?
%  ! check quickly whether superimposed trials at FCz look better with ASR5 or 10 and also check averages/t-tests (afterwards!)
%  ! add ICA computation after merging blocks
%     ! make sure to plot topos properly when evaluating ICs
%     !  need to remove non-EEG channels first
%     ! since I have chanlocs loaded in I may as well try to run the ICLabel algorithm to see how certain it is about potential eye/muscle artifacts
%
%
%
%%
clc
clearvars
% addlab

% subfold = 'main'; % main experiment data
subfold = 'pilots_techtests';

homedir = {'/home/rick/neuro/iannettilab/Projects/VP monkey/data/raw_' subfold];
cd(homedir);

addpath('/home/rick/neuro/iannettilab/Projects/VP monkey/Giac_ToolBox')

%% SETTINGS
s = [];
s.savePath =  ['/home/rick/neuro/iannettilab/Projects/VP monkey/data/clean_' subfold]; % mkdir(s.savePath)

% conditions
s.conds = {'VIS', 'SOM', 'AUD'};

% re-referencing
s.chans4ref   = {'EXG1','EXG2'}; s.refname = 'ears'; % EXG channels to use as reference for all EEG channels
% s.chans4ref = [1:28]; s.refname = 'CAR'; % avg ref
% s.chans4ref = {'AF4','FC6','PO4','POz','PO3','FC5','AF3'}; s.refname = 'ring'; 
s.chans2excludefromref = {'Erg1'};
s.chans2remove  = { 'EXG3','EXG4','EXG5','EXG6','EXG7','EXG8', 'X1','X2','X3','X4', 'Status'   };
s.chans2replace = {'Cz','CPz', 'Fz', 'F2' }; % missing electrodes in monkey caps

% ASR
s.asr.useASR = true;
s.asr_useRiemmanian = false; % faster I think
s.asr_useGPU  = false;
s.asr_thresh  = 5;
s.asr_MaxMem  = 16 * 1024; % MB, using automatic now

% frequency filter
s.freqs       = [1, 30]; % band pass frequencies

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi

%% select sessions
sesh = dir; sesh = {sesh.name}; sesh(1:2) = [];
sel = listdlg('ListString', sesh, 'SelectionMode','multiple');
% sel = 2; warning( '! BYPASSING SESSION SELECTION FOR SESSION P02' )
sesh = sesh(sel);
nsesh = length(sesh);

%% loop through sessions
for sh = 1 :nsesh

    % grab raw data file names
    cd([homedir filesep sesh{sh}])
    files = dir; files = {files.name}; files(~endsWith(files, '.bdf')) = [];
    nfiles = length(files);

    %% loop through files (blocks)
    tIN_file = tic;
    EEG_all = []; % all for this session
    numChansInterpolated = 0;  % counting for later ICA
    for f = 1:nfiles

        % load file 
        if      endsWith( files{f}, '.bdf' )
%           EEG = pop_biosig( files{f} , 'importevent','off', 'rmeventchan' ,'off')
            EEG = pop_biosig( files{f} , 'importevent','off', 'rmeventchan' ,'off');
        elseif  endsWith(files{f}, '.set')
            EEG = pop_loadset(files{f}); % if raw continuous data have already been imported / processed 
        end

        % extract triggers
        trigChan = EEG.data(strcmp({EEG.chanlocs.labels}, 'Status'),:);
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
%         if length(triggers) > 0
            events = struct;
            for trl = 1 : length(triggers)
                events(trl).type = s.conds{trigCodes(trl)};
                events(trl).urevent =  trl;
                events(trl).latency =  triggers(trl);
            end
%         else
%             events = h.events;
%         end
        EEG.event = events; clear events % add events to EEG structure
        EEG.urevent = EEG.event; EEG.urevent = rmfield(EEG.urevent, 'urevent');
        EEG = eeg_checkset(EEG); 
        
        % remove channels we don't need
        EEG = pop_select( EEG, 'nochannel',  find( cellfun( @(x) ismember(x, s.chans2remove) , {EEG.chanlocs.labels} ) ) );

        % select only EEG channels for main preprocessing steps
        chans_EEG_EXG = strcmp({EEG.chanlocs.type}, 'EEG' ) | startsWith({EEG.chanlocs.labels}, 'EXG' )  ;
        EEG_otherchans  = pop_select( EEG, 'channel',  find(~chans_EEG_EXG));
        EEG             = pop_select( EEG, 'channel',  find( chans_EEG_EXG));
        
        % band-pass filter EEG & EXG channels  
        EEG.data = EEG.data - mean(EEG.data,2);  % remove DC offset
        EEG = pop_eegfiltnew( EEG, 'locutoff', s.freqs(1)); % high-pass
        EEG = pop_eegfiltnew( EEG, 'hicutoff', s.freqs(2)); % low-pass
%         EEG = pop_firws(EEG, 'wtype','hamming',  'fcutoff', s.freqs,  'forder', pop_firwsord('hamming', EEG.srate, 0.5) ); % 'usefftfilt',true, 'minphase',true 

        % remove EXG chans from EEG and add to EEG_otherchans (so they are excluded from ASR)
        EEG_EXGchans = pop_select( EEG, 'channel', find( startsWith({EEG.chanlocs.labels}, 'EXG' )) );
        EEG          = pop_select( EEG, 'channel', find(~startsWith({EEG.chanlocs.labels}, 'EXG' )) );
        EEG_otherchans = rs_eeglab_mergechans(EEG_EXGchans, EEG_otherchans);

        % add chanlocs to EEG chans
        for c = 1:EEG.nbchan
            temp = chanlocs(c);
            if ~strcmp( temp.labels, EEG.chanlocs(c).labels ), error 'channel labels do not match!!'; end
            temp = rmfield(temp,'labels');
            flds = fields(temp);
            for k = 1:length(flds)
                EEG.chanlocs(c).(flds{k}) = temp.(flds{k});
            end
        end
        
        % clean with ASR --------
        if s.asr.useASR
            fprintf('running ASR cleaning ...\n')
            % find clean data for calibration
            [ref_section, ref_mask] = clean_windows(EEG); % , s.ref_maxbadchannels, s.ref_tolerances, s.ref_wndlen);
            % apply ASR cleaning 
            tin_asr = tic;
            EEG2clean = EEG; 
            rng('default');
            EEG2clean = clean_asr( EEG2clean, s.asr_thresh, [],[],[],     ...
                           ref_section,[],[], s.asr_useGPU , s.asr_useRiemmanian, s.asr_MaxMem ); 
            s.asr_tElapsed = toc(tin_asr);
            fprintf( '... cleaning completed in %.2f s\n\n', s.asr_tElapsed ) 
    %         vis_artifacts(EEG2clean, pop_select(EEG,'channel',1:64)); % visualise data before and after cleaning
            EEG = EEG2clean; clear EEG2clean
        else
            disp 'skipping ASR ...'
        end

        % epoch data
        EEG             = pop_epoch(EEG, s.conds, [-2 2 ]);
        EEG_otherchans  = pop_epoch(EEG_otherchans, s.conds, [-2 2 ]);

        % interpolate missing channels (s.chans2replace) & bad channels from their neighbours
        % find bad channels using Giac function
        dataft = eeglab2fieldtrip( EEG, 'raw', 'none' ); 
        badchans = Giac_EEG_CatchNoisyElectrodes(dataft, {'all'}, 3, 'recursive');
        chans2interpolate = unique([badchans',  s.chans2replace]);
        chans2interpolate = cell2mat( cellfun(@(x) find(strcmp({EEG.chanlocs.labels}, x)), chans2interpolate, 'UniformOutput', false) );
        numChansInterpolated = max( [length(chans2interpolate), numChansInterpolated]);

        %!!! interpolate within fieldtrip
        cfg = []; cfg.elec = dataft.elec; lay = ft_prepare_layout(cfg); clear cfg; % prepare layout using the chanlocs gained from locs file

        EEG = pop_interp(EEG, chans2interpolate, 'spherical'); % ? compare/investigate this sexy 'spacetime' option which interpolates both in space and time
                                                               % !  compare also to Giacomo's approach. e.g. this method doesn't ask for neighbours ...
        % add other chans back in 
        EEG = rs_eeglab_mergechans(EEG, EEG_otherchans);

        % re-reference all EEG channels        
        if isnumeric(s.chans4ref)
            s.chans4ref = {EEG.chanlocs(s.chans4ref).labels}; clear labs
        end
        EEG = pop_reref(EEG, s.chans4ref, 'keepref', 'on', ...
            'exclude' , find(cellfun(@(x)ismember(x,s.chans2excludefromref),{EEG.chanlocs.labels}))); 
        
        % store data for later merging
        EEG_all = [EEG_all,  EEG];

        fprintf('finished file %d/%d (tElapsed = %.2f mins)\n', f, nfiles, toc(tIN_file)/60  )
    end % end file (block) loop

    % merge epoched files
    EEG = pop_mergeset( EEG_all, 1:nfiles ); 
    
    % compute ICA
     % !remove non-EEG electrodes? or keep EXG1 and EXG2 
     % !num IC should be numchans - numChansInterpolated

    % save
    if s.asr.useASR, ASRname = num2str(s.asr_thresh); else, ASRname  = 'no'; end
    saveName = [ 'ASR_' ASRname ' reref_' s.refname ];
    savePath = [s.savePath filesep sesh{sh} filesep]; mkdir(savePath)
    filename = [ saveName  ' merged ' sesh{sh}];
    EEG.setname = filename;
    pop_saveset(EEG, 'filepath', savePath, 'filename', filename  );
    fprintf('... saved session %s: %s\n', sesh{sh}, [savePath filename]) 

    % export to letswave (preserve ICA weights if possible!)
    %!!




end % end session loop

