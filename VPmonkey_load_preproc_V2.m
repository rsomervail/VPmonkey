%   
% 
%   
%%
clc
clearvars
close all

subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';

homedir = ['/media/rick/Rick_LabData_3/neuro/iannettilab_ongoing/VPmonkey/data/raw_' subfold];
cd(homedir);

suadir = '/media/rick/Rick_LabData_3/neuro/iannettilab_ongoing/VPmonkey/data/raw_SUA';

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

% try 
%     eeglab
% catch
%     addlab
% end

%% SETTINGS
s = [];
saveas = '';
% saveas = ['_SPK_thresh_45']
s.savePath = [ strrep(homedir,'raw','clean') saveas]; 
if ~exist( s.savePath,'dir'), mkdir(s.savePath); end
if ~exist( [s.savePath filesep 'lw'],'dir'), mkdir([s.savePath filesep 'lw']); end

% choose modalities to process 
s.eeg  = false;
s.lfp  = false;
% s.mua  = false; % OLD method computing MUA from electrodes myself
s.sua  = true; % handles MUA as well
s.eye  = false;
s.misc = false;

% session list
s.seshlist = [getRoot '/VPmonkey/VPmonkey_sessionList.csv'];
seshlist = VPmonkey_importSessionlist; %(s.seshlist, [4, Inf]);

% conditions
s.conds = {'VIS', 'SOM', 'AUD'};

% baseline correction
s.blwin_eeg = [  -200   -20  ]; % EEG / LFP
s.blwin_pup = [  -500   -20  ]; % pupil diameter (probs would need a shorter one for gaze direction if I need that at some point)
s.blwin_spk = [  -1500  -500  ]; % MUA & SUA

% MUA & SUA kernel standard deviation in ms
s.sdf_kernelSD = 30;  % ? crazy small but need to try an extreme value just in case

% use data visualisation
s.plotData = false;
s.plotYlims = [-100 100]; % plot limits for ft_databrowser (interpolating chans)
s.plotYlims2 = [-15 15]; % plot limits for ft_databrowser (rejecting epochs)

% epoching 
s.xlims.preproc  = [-2 8]; % initial limits for all data
s.xlims.lab.EEG  = [-0.5 2];  
s.xlims.lab.LFP  = [-0.5 2];  
% s.xlims.lab.MUA  = [-1 2];   % this is determined by the epoched range of spikes I get
s.xlims.lab.MISC = [-0.5 0.5];
s.xlims.lab.EYE  = [-0.5 8]; 
s.xlims.lw.EEG   = [-0.5 1];  
s.xlims.lw.LFP   = [-0.5 1];  
s.xlims.lw.MUA   = [-1 1];  
s.xlims.lw.MISC  = [-0.2 0.5];
s.xlims.lw.EYE   = [-0.5 6]; 

% downsampling factor for lw export
s.dsf = 2; 

% re-referencing
% s.chans4ref   = {'EXG1','EXG2'}; s.refname = 'ears'; % EXG channels to use as reference for all EEG channels
% s.chans4ref = [1:28]; s.refname = 'CAR'; % avg ref
s.chans4ref = {'AF4','FC6','PO4','POz','PO3','FC5','AF3'}; s.refname = 'ring'; 
s.chans2remove  = { 'EXG4','EXG5','EXG6','EXG7','EXG8', 'X1','X2','X3','X4', ...
        'Status', 'GSR1', 'GSR2', 'Erg2', 'Resp', 'Plet', 'Temp'  };

% ASR
s.asr.useASR = true;
s.asr_useRiemmanian = false;  % Riemmannian faster but not independently validated
s.asr_useGPU  = false; %? not necessarily faster, especially if can use lots of RAM for CPU
s.asr_thresh  = 15;
s.asr_MaxMem  = 24 * 1024; % MB, using automatic now

% frequency filter
s.freqs       = [1, 98]; % band pass frequencies

s.lfp_conds = {'VIS', 'SOM', 'AUD'};
s.lfp_trig      = [2, 8, 32];

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % REORDERED BECAUSE CHANLOCS REVERSED IN THIS FILE
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw); % ? once had issues because I changed the variable name in this file to chanlocs and overwrote the other one

% interpolation
s.int_neighdist = 0.25;  % for interpolation 

%% import intracranial information table
tbl_depths = VPmonkey_import_electrodeDepths;

%% select sessions
sesh = dir; sesh = {sesh.name}; sesh(1:2) = [];
sel = listdlg('ListString', sesh, 'SelectionMode','multiple');
% sel = 2; warning( '! BYPASSING SESSION SELECTION' )
sesh = sesh(sel);
nsesh = length(sesh);

%% loop through sessions
for sh = 1 :nsesh

    % get sesh / sub info
    subid   = sesh{sh}(end);
    seshnum = str2double(sesh{sh}(2:3));

    % get subject-specific missing electrodes (currently the same, sometimes FCz is also bad in either mk)
    if      strcmp( subid ,'M')
        s.chans2replace = {'CZ','CPz', 'FCz', 'F2' }; 
    elseif  strcmp( subid, 'T')
        s.chans2replace = {'CZ','CPz', 'Fz',  'F2' };
    end

    % get session-specific modality offsets
    seshrow = find(strcmp( seshlist.session ,sesh{sh}(1:3)));
    s.cond_offsets = [ seshlist.vis_offset(seshrow) , seshlist.som_offset(seshrow), seshlist.aud_offset(seshrow)    ]; % ms

    % grab raw data file names
    cd([homedir filesep sesh{sh}])
    files = dir; files = {files.name}; 
    files(~endsWith(files, '.bdf') & ~endsWith(files, '.set')) = [];
    files(startsWith(files,'SPK')) = [];
    nfiles = length(files);

    %% loop through files (blocks)
    tIN_file = tic;
    EEG_all = []; % all for this session
    numChansInterpolated = 0;  % counting for later ICA
    for f = 1 :nfiles

        %% load file 
        if      endsWith( files{f}, '.bdf' )
            try 
                EEG = pop_biosig( files{f} , 'importevent','off', 'rmeventchan' ,'off', 'channels', 1:48); % w/  EXG
            catch
                EEG = pop_biosig( files{f} , 'importevent','off', 'rmeventchan' ,'off', 'channels', 1:40); % w/o EXG
            end
        elseif  endsWith(files{f}, '.set')
            EEG = pop_loadset(files{f}); % if raw continuous data have already been imported / processed 
%             if s.mua
%                 SPK = pop_loadset(['SPK ' files{f}]);
%             end
        end

        %% file-wise preprocessing steps 
        
        % get events if not already done in a previous run of the dual EEG splitting script
        if isempty(EEG.event)
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
        end

%         % if microphone is present, rectify the signal
%         if any(strcmp({EEG.chanlocs.labels},'micro'))
%             chanMicro = find(strcmp({EEG.chanlocs.labels},'micro'));
%             EEG.data( chanMicro,:) = abs(EEG.data( chanMicro,:)); % rectify
%         end

        %% OLD - MUA processing (before new spikes from spike sorting analysis)
%         if s.mua 
%             fprintf('processing spikes ...\n')
% 
%             % scale units of spike data to microvolts (to match EEG)    
%             SPK.data = SPK.data * 1e6; % ? not strictly necessary but nicer to know everything's the same
%     
%             % band-pass filter  ! should ideally replace with a wavelet filter .. may need to code from scratch using Mike Cohen book
% %             SPK = pop_eegfiltnew( SPK, 'locutoff',  300); % high-pass
%             SPK = pop_eegfiltnew( SPK, 'locutoff',  300, 'hicutoff', 5000); % band-pass (recommended by Eros)
% 
%             % define MUA structure 
%             MUA = pop_select(EEG,'channel',find(strcmp({EEG.chanlocs.type},'LFP')));
%             for k = 1:5
%                 MUA.chanlocs(k).labels = strrep( MUA.chanlocs(k).labels, 'LFP','MUA' );
%                 MUA.chanlocs(k).type = 'MUA';
%             end
%             MUA.data = zeros(size(MUA.data));
% 
%             % mark bad channels with NaNs
%             depth = tbl_depths.Depth(strcmp(tbl_depths.sub,subid) ...
%                 & strcmp(tbl_depths.sesh,['S' sprintf('%02d',seshnum)]));
%             MUA.data(isnan(depth),:) = nan;
%       
%             % threshold spikes  ! have maximum amplitude too? since there may be some high-frequency artifacts ... although maybe they can just be really close. FIND OUT APPROPRIATE RANGE
%             thresh = 45;
%             minspikedist = 1/1000; 
%             for k = 1:SPK.nbchan
% 
%                 % if bad channel then may as well skip this step
%                 if ~any(isnan(MUA.data(k,:)))
% 
%                     % get latencies of putative spikes
%                     [pks, locs] = findpeaks( -SPK.data(k,:), "MinPeakHeight",thresh,...
%                             "MinPeakDistance", round(minspikedist * SPK.srate) ...
%                             );
%                     lats = SPK.times(round(locs))/1000; 
%     
%                     % if any putative spikes were found...
%                     if ~isempty(lats)
%     
%                         % compute spike density function
%                         sdf = f_calcSDF(lats, s.sdf_kernelSD, EEG.times(end)/1000 );
%                         sdf(end) = [];
%                         sdf_times = 0 : (1/1000) : (EEG.times(end)/1000);
%     
%                         % resample to match EEG
%                         sdf_rs = interp1(sdf_times, sdf, EEG.times/1000);
%     
%                         MUA.data(k,:) = sdf_rs;
%     
%                     end % end is empty lats
% 
%                 end % end check for bad channel 
%                
%             end % end loop through chans
% 
%             % ! would be superceded by having the true spikes from spike-sorting
% %             % linear interpolation to remove shock artifact  ?? might need low-pass filtering, or longer movmean window for this
% %             ctemp = [];
% %             ctemp.eventcodes  = {'SOM'};
% %             ctemp.timewindow = [ 10 10 ]/1000; 
% %             ctemp.channels =  1:MUA.nbchan;
% %             MUA = pop_rs_interpolateShockArtifacts(MUA,ctemp);
%             
%             % epoch
%             MUA = pop_epoch(MUA, s.conds, s.xlims.preproc);
% 
%             % compute mean of MUA channels      
%             % ? this probably reflects population code more, but ignores any differences in depth
%             MUA.data(end+1,:,:) = mean(MUA.data,'omitnan');
%             MUA.chanlocs(end+1) = MUA.chanlocs(end);
%             MUA.chanlocs(end).labels = 'MUA_avg';
%             MUA.nbchan = MUA.nbchan + 1;
% 
%             fprintf('... finished processing spikes\n')
%         end
%         % END SPIKE PROCESSING 

        %% LFP - scale units of LFP to microvolts (to match EEG)
        if s.lfp
            chans_lfp = strcmp({EEG.chanlocs.type},'LFP');
            EEG.data(chans_lfp,:) = EEG.data(chans_lfp,:) * 1e6;
        end

        %% notch filter LFP
        if s.lfp
            chans_lfp = find(strcmp({EEG.chanlocs.type},'LFP'));
            temp = pop_select(EEG, 'channel', chans_lfp);
            temp = pop_eegfiltnew( temp, 'locutoff', 48, 'hicutoff', 52, 'revfilt', true); % notch 
            % ? add harmonics at 100, 150, 200 & 250?
            EEG.data(chans_lfp,:) = temp.data;
        end

        %% LFP - mark bad channels with NaNs
        if s.lfp
            % mark bad channels with NaNs
            depth = tbl_depths.Depth(strcmp(tbl_depths.sub,subid) ...
                & strcmp(tbl_depths.sesh,['S' sprintf('%02d',seshnum)]));
            for k = 1:5
                if isnan(depth(k))
                    channum = find(strcmp({EEG.chanlocs.labels},['LFP_el' num2str(k)]));
                    EEG.data(channum,:) = nan;
                end
            end
        end

        %% LFP - take average of channels
        if s.lfp
            % compute mean of LFP channels      
            % ? this probably reflects population activity more, but ignores any differences in depth
            lfpchans = find(strcmp({EEG.chanlocs.type},'LFP'));
            EEG.data(end+1,:)   = mean(EEG.data(lfpchans,:),'omitnan');
            EEG.chanlocs(end+1) = EEG.chanlocs(end);
            EEG.chanlocs(end).labels = 'LFP_avg';
            EEG.nbchan = EEG.nbchan + 1;
        end

        %% correct latencies of stimulus events if mean error is specified        ? I left some comment about issue with importing latency error from spreadsheet but seemed fine to me now ...
        offsetSamps = round((s.cond_offsets/1000) * EEG.srate); % assumes offsets are in ms here
        for trl = 1:length(EEG.event)
            if offsetSamps(strcmp(EEG.event(trl).type,s.conds))~=0
                EEG.event(trl).latency = EEG.event(trl).latency + offsetSamps(strcmp(EEG.event(trl).type,s.conds));
            end
        end
              
        %% linearly interpolate any shock artifacts in selected channels
        if s.eeg
            ctemp = [];
            ctemp.eventcodes  = {'SOM'};
            ctemp.timewindow = [ 5 10 ]/1000; 
            ctemp.channels =  find(strcmp({EEG.chanlocs.type},'EEG'));
    %         ctemp.channels =  find(strcmp({EEG.chanlocs.type},'EEG') | strcmp({EEG.chanlocs.type},'LFP'));
            EEG = pop_rs_interpolateShockArtifacts(EEG,ctemp);
        end

        %% linearly interpolate any visual stimulus artifacts in pupil channels
        if s.eye
            ctemp = [];
            ctemp.eventcodes  = {'VIS'};
            ctemp.timewindow = [ 100 150 ]/1000; 
            ctemp.channels =  find(startsWith({EEG.chanlocs.labels},'eye_pupil','ignorecase',true));
            EEG = pop_rs_interpolateShockArtifacts(EEG,ctemp);
        end

        %% pupil - Savitsky-Golay filter (smoothing & removal of dropouts)
%         if s.eye
% 
%             % get pupil values
%             pupW = double(EEG.data( strcmp({EEG.chanlocs.labels},'eye_PupilWidth'),  :)); %   figure; plot(pupW)
%             pupH = double(EEG.data( strcmp({EEG.chanlocs.labels},'eye_PupilHeight'), :)); %   figure; plot(pupH)
%             pupD = mean( [pupW; pupH] ); % raw pupil diameter
%     
%     %         % exclude extreme values
%     %         thresh = 0.04; 
%     %         pupW(pupW<thresh) = nan;  
%     %         pupH(pupH<thresh) = nan; 
% 
%             % savitsky-golay filter  %? what parameters to choose here
%             order    = 1;
%             framelen = floor(0.5 * EEG.srate);
%             if mod(framelen,2) == 0
%                 framelen = framelen-1;
%             end
%             pupW_filt = sgolayfilt(pupW,order,framelen);
%             pupH_filt = sgolayfilt(pupH,order,framelen);
%             pupD_filt = mean( [pupW_filt; pupH_filt]  ); % DON'T omitnan here or there are sudden jumps caused by one channel popping in suddenly after a nan period
%             
%             % plot
%             t = EEG.times/1000;
%             figure; subplot(2,1,1); plot(t,pupW); hold on; plot(t,pupW_filt);
%                     subplot(2,1,2); plot(t,pupH); hold on; plot(t,pupH_filt);
% %             figure; plot(t,pupD); hold on; plot(t,pupD_filt);
% 
% %             % Replace old pupil diameter with preprocessed version
% %             chan = find(strcmp({EEG.chanlocs.labels},'eye_pupilDiameter'));
% %             EEG.data( chan, :) = pupD_filt;
% %             EEG.chanlocs(chan).labels = 'eye_PupilDiameter'; % rename for consistency
% %     
% %             % Replace pupil height with raw pupil diameter
% %             chan = find(strcmp({EEG.chanlocs.labels},'eye_PupilHeight'));
% %             EEG.data( chan, :) = pupD;
% %             EEG.chanlocs(chan).labels = 'eye_PupilDiameterRAW';
% 
% 
% 
%         end


        %% pupil - moving median low-pass filter (smoothing & removal of dropouts)
        if s.eye

            % get pupil values
            pupW = EEG.data( strcmp({EEG.chanlocs.labels},'eye_PupilWidth'),  :); %   figure; plot(pupW)
            pupH = EEG.data( strcmp({EEG.chanlocs.labels},'eye_PupilHeight'), :); %   figure; plot(pupH)
            pupD = mean( [pupW; pupH] ); % raw pupil diameter
    
    %         % exclude extreme values
    %         thresh = 0.04; 
    %         pupW(pupW<thresh) = nan;  
    %         pupH(pupH<thresh) = nan; 
    
            % compute moving median
            win = 0.20; % seconds (? should not be greater than 1 second as main pupil component seems to be in the 0-1Hz range)
            pupW_filt = movmedian( pupW , win*EEG.srate, 'omitnan' );  %   figure; plot(pupW); hold on; plot(pupW_filt)
            pupH_filt = movmedian( pupH , win*EEG.srate, 'omitnan' );  %   figure; plot(pupH); hold on; plot(pupH_filt)
            pupD_filt = mean( [pupW_filt; pupH_filt]  ); % DON'T omitnan here or there are sudden jumps caused by one channel popping in suddenly after a nan period
%             figure; plot(pupD); hold on; plot( pupD_filt );

            % Replace old pupil diameter with preprocessed version
            chan = find(strcmp({EEG.chanlocs.labels},'eye_pupilDiameter'));
            EEG.data( chan, :) = pupD_filt;
            EEG.chanlocs(chan).labels = 'eye_PupilDiameter'; % rename for consistency
    
            % Replace pupil height with raw pupil diameter
            chan = find(strcmp({EEG.chanlocs.labels},'eye_PupilHeight'));
            EEG.data( chan, :) = pupD;
            EEG.chanlocs(chan).labels = 'eye_PupilDiameterRAW';
        
        end
        
        %% split EEG/EXG chans from others
        chans_EEG_EXG = strcmp({EEG.chanlocs.type}, 'EEG' ) | startsWith({EEG.chanlocs.labels}, 'EXG' )  ;
        EEG_otherchans  = pop_select( EEG, 'channel',  find(~chans_EEG_EXG));
        EEG             = pop_select( EEG, 'channel',  find( chans_EEG_EXG));
        
        %% band-pass & notch filter EEG & EXG channels  
        if s.eeg
            EEG.data = EEG.data - mean(EEG.data,2);  % remove DC offset
            EEG = pop_eegfiltnew( EEG, 'locutoff', s.freqs(1)); % high-pass
            EEG = pop_eegfiltnew( EEG, 'hicutoff', s.freqs(2)); % low-pass
            EEG = pop_eegfiltnew( EEG, 'locutoff', 48, 'hicutoff', 52, 'revfilt', true); % notch 
    %         EEG = pop_firws(EEG, 'wtype','hamming',  'fcutoff', s.freqs,  'forder', pop_firwsord('hamming', EEG.srate, 0.5) ); % 'usefftfilt',true, 'minphase',true 
        end

        %% find blinks using EEG and pupil data and store in event structure
        if s.eye 

            % if not doing most of EEG preprocessing, need to quickly do the high-pass filter at least
            if ~s.eeg 
                EEG.data = EEG.data - mean(EEG.data,2);  % remove DC offset
                EEG = pop_eegfiltnew( EEG, 'locutoff', s.freqs(1)); % high-pass
            end
            
            % find blinks in EEG & pupil        
            % PUPIL
            pup = -EEG_otherchans.data(strcmp({EEG_otherchans.chanlocs.labels},'eye_PupilDiameterRAW'),:);
            thresh = -0.06;
            minpeakdistance = 0.2;
            minwidth = 20/1000;
    
            [~,lats, W] = findpeaks(pup,'MinPeakDistance', minpeakdistance*EEG.srate , 'MinPeakHeight', thresh  ...
                ,'MinPeakWidth', round(minwidth*EEG.srate) );
    %             ,'MaxPeakWidth', round(maxwidth*EEG.srate) );
            widths = W / EEG.srate; clear W
    
            % EEG
            thresh_eeg = 2; % z-score
            minwidth_eeg = 20 /1000;
            %
            chans = {EEG.chanlocs.labels};
            chaninds = find( (startsWith(chans,'F') | startsWith(chans,'AF')) & ~startsWith(chans,'FC')  &  ~ismember(chans,s.chans2replace));
            eeg = EEG.data( chaninds, :);
            eeg = movmean(eeg, round( EEG.srate * 1/15), 2  ); % quick low-pass filter to make blink identification more robust
            eeg = mean(eeg); % ? taking mean of all these channels may not be most sensitive approach
            eeg = (eeg - mean(eeg)) / std(eeg);  
            [~,lats_eeg, W] = findpeaks(eeg,'MinPeakDistance', minpeakdistance*EEG.srate , ...
                'MinPeakHeight', thresh_eeg ,'MinPeakWidth', round(minwidth_eeg*EEG.srate) );
            widths_eeg = W / EEG.srate; clear W  
    
            % combine EEG- and pupil- identified blinks
            mindist = 0.1 * EEG.srate; % EEG and pupil peaks shouldn't be too far from each other
            inds2remove = [];
            for k = 1:length(lats)
                if ~any( abs(lats(k) - lats_eeg) < mindist )
                    inds2remove = [inds2remove, k];
                end
            end
            lats(inds2remove) = [];
    
            % linearly interpolate over all identified blinks
            if ~isempty(lats)
                EEG_otherchans = eeg_addnewevents(EEG_otherchans, num2cell(lats), repmat({'BLINK'},size(lats)));
                
                ctemp = [];
                ctemp.eventcodes  = {'BLINK'};
                ctemp.timewindow = [170; 300]/1000; 
                ctemp.timewindow = [170; 300]/1000; 
                ctemp.channels = find(strcmp({EEG_otherchans.chanlocs.labels},'eye_PupilDiameter'));
                EEG_otherchans = pop_rs_interpolateShockArtifacts(EEG_otherchans,ctemp);
                EEG_otherchans.event( strcmp({EEG_otherchans.event.type},'BLINK') ) = [];
            end

%         pop_letsplot([temp,EEG])
%         clear temp lats pup

        end


        %% handle channels & channel locations
        % if any are present, remove EXG chans from EEG and add to EEG_otherchans (so they are excluded from ASR)
        if any(startsWith({EEG.chanlocs.labels}, 'EXG' ))
            EEG_EXGchans = pop_select( EEG, 'channel', find( startsWith({EEG.chanlocs.labels}, 'EXG' ) ) );
            EEG          = pop_select( EEG, 'channel', find(~startsWith({EEG.chanlocs.labels}, 'EXG' ) ) );
            EEG_otherchans = rs_eeglab_mergechans(EEG_EXGchans, EEG_otherchans); clear EEG_EXGchans
        end

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
        
        %% clean with ASR -------- 
        if s.asr.useASR && s.eeg
                fprintf('running ASR cleaning ...\n')
            % make copy of EEG data without the missing channels (but not removing the noisy channels to  be interpolated which are found after ASR cleaning)
                tin_asr = tic;
            EEG2clean = pop_select(EEG, 'nochannel', s.chans2replace); % remove missing channels from ASR computation for better sensitivity
            % find clean data for calibration
            [ref_section, ref_mask] = clean_windows(EEG2clean); % , s.ref_maxbadchannels, s.ref_tolerances, s.ref_wndlen);
            % check enough data were found for good calibration
            if      ref_section.times(end)/1000 < 15
                error('NOT ENOUGH DATA FOR ASR CALIBRATION IN THIS BLOCK (%.1fs) - CONSIDER REJECTING?', ref_section.times(end)/1000)
            elseif  ref_section.times(end)/1000 < 60
                warning('SMALL AMOUNT OF DATA FOR ASR CALIBRATION (%.1fs), BUT PROBABLY OK',ref_section.times(end)/1000)
            end
                rng('default'); % can't remember why this was here 
            % apply ASR cleaning 
            EEGcleaned = clean_asr( EEG2clean, s.asr_thresh, [],[],[],     ...
                           ref_section,[],[], s.asr_useGPU , s.asr_useRiemmanian, s.asr_MaxMem ); 
                s.asr_tElapsed = toc(tin_asr);
                fprintf( '... cleaning completed in %.2f s\n\n', s.asr_tElapsed ) 
                    if s.plotData
                        vis_artifacts(EEGcleaned, EEG2clean, 'ScaleBy', 'new' ); % visualise data before and after cleaning  % !!! NEED A WAY TO CLOSE THIS AFTER CHECKING DATA, SO THAT IT DOESN'T FUCK WITH GIAC INTERP PLOTTING
                        input('press enter when ready to continue\n','s'); % !!! improve this by some kind of wait for figure close (how does ft_databrowser do it?)
                        close all
                    end
            % put cleaned data channels back into main EEG structure
            for c = 1:EEGcleaned.nbchan
                chanindy = find(strcmp({EEG.chanlocs.labels},EEGcleaned.chanlocs(c).labels));
                EEG.data(chanindy,:) = EEGcleaned.data(c,:);
            end
            clear EEG2clean EEGcleaned ref_section
        else
            disp 'skipping ASR ...'
        end

        %% epoch data
        EEG             = pop_epoch(EEG,            s.conds, s.xlims.preproc);
        EEG_otherchans  = pop_epoch(EEG_otherchans, s.conds, s.xlims.preproc);

        %% interpolate bad channels from their neighbours -----------------------
        if s.eeg
            % convert to fieldtrip
            EEGtemp = pop_select(EEG, 'nochannel', find( cellfun(@(x) ismember(x,s.chans2replace), {EEG.chanlocs.labels}) | ~strcmp({EEG.chanlocs.type},'EEG')) );  % remove missing channels & non-EEG channels from fieldtrip data
            dataft = eeglab2fieldtrip( EEGtemp, 'raw', 'none' ); clear EEGtemp  % convert to fieldtrip
            for k = 1:length(dataft.trial), dataft.trial{k} = double(dataft.trial{k} ); end % interpolation functions used below require double-precision data
            % interpolate bad channels (within fieldtrip)
            [ dataft, interpolated_channels ] = RS_Giac_EEGinterpolation( dataft, 'all', s.plotYlims, s.chanlay, s.int_neighdist, s.plotData );
    
            % interpolate missing channels (within fieldtrip)
            dataft = RS_Giac_EEG_CreateElectrodes(dataft, strrep( s.chans2replace, 'CZ', 'Cz' ), s.chanlay, .20); % .20 is better for central electrodes
            numChansInterpolated = max( [length(interpolated_channels)+length(s.chans2replace), numChansInterpolated]);
    
            % convert data back into EEGlab (simply replace EEG.data with data from interpolated fieldtrip version)
            EEG_interped = fieldtrip2eeglab( dataft ); chansft = dataft.label; clear dataft
            for c = 1:EEG_interped.nbchan
                chanindy = find(strcmpi( chansft{c}, {EEG.chanlocs.labels} )); % get the correct index for this channel (order is different in fieldtrip data)
                EEG.data(chanindy,:,:) = EEG_interped.data(c,:,:); % replace original EEG data with data after interpolation in fieldtrip
            end; clear chansft EEG_interped chanindy
        end
        
        %% add other chans (i.e. EXG, Erg1, eyetracker or LFP channels) back in 
        EEG = rs_eeglab_mergechans(EEG, EEG_otherchans); clear EEG_otherchans

        %% re-reference all EEG channels        
        if s.eeg
            if isnumeric(s.chans4ref), s.chans4ref = {EEG.chanlocs(s.chans4ref).labels}; end % ? ref channels needs to be cell array
            EEG = pop_reref(EEG, s.chans4ref, 'keepref', 'on', ...
                 'exclude', find( ~strcmp({EEG.chanlocs.type},'EEG') ) ); % exclude non-EEG channels from re-reference 
    %             'exclude' , find(cellfun(@(x)ismember(x,s.chans2excludefromref),{EEG.chanlocs.labels}))); 
        end

        %% OLD - append MUA channels for easier merging (before new SUA integration)
%         if s.mua
%             EEG = rs_eeglab_mergechans(EEG, MUA);
%         end
        
        %% store data for later merging
        EEG_all = [EEG_all,  EEG]; 

        fprintf('finished file %d/%d (tElapsed = %.2f mins)\n\n\n', f, nfiles, toc(tIN_file)/60  )
    end % end file (block) loop

    %% merge epoched files
    if length(EEG_all) > 1
        EEG = pop_mergeset( EEG_all, 1:nfiles );
        EEG.urevent( contains({EEG.urevent.type},'boundary') ) = []; % remove boundary events to make trial rejection stuff easier later
        for k = 1:length(EEG.event)
            EEG.event(k).urevent = k; % then match each event to the new set of urevents
        end
    else
        EEG = EEG_all; % if only one file in session then no need to merge
    end; clear EEG_all

    %% remove any repeat events (due to epochs being long enough to overlap with other stimulus events)
    EEG = pop_rs_removeNonZeroEvents(EEG); 

    %% LOAD SUA DATA AND INTEGRATE
    % ? Trials.SUA 1st column is seconds, absolute time, 2nd col is electrode number (1:5 = SubT, 6:10 = SubM), 3rd col is unit
    %  unit nums may be repeated for different electrodes, so should store electrode num in SUA channel labels too 
    % ? can compare spike times to stimulus time by looking at Trials.event, 
    %  2nd event is the stimulus, code is 2, 8, 32 (Vis, Som, Aud)
    if s.sua        
        cd(suadir)

        % find SUA for this subject
        file2load = ['S' sprintf('%02d',seshnum) '_MT_SUA.mat' ];
        load(file2load);
        sua = Trials; clear Trials

        % check number of events match
        if EEG.trials ~= length(sua)
            error 'TRIAL NUMBER MIS-MATCH BETWEEN SUA AND EEG'
        end

        % take only the electrodes for this subject
        if strcmp(subid,'T')
            elec2take = 1:5;
        elseif strcmp(subid,'M')
            elec2take = 6:10;
        end
        units = [];
        for trl = 1:length(sua)
            temp = sua(trl).SUA;
            inds = ismember(temp(:,2),elec2take);
            sua(trl).SUA = temp(inds,:);
            % for subject M need to subtract 5 from electrode number
            if strcmp(subid,'M')
                sua(trl).SUA(:,2) = sua(trl).SUA(:,2) - 5;
            end
            units = [units; sua(trl).SUA(:,2:3)];
        end


        % get complete list of units
        units2 = cell(length(units),1);
        for k = 1:length(units)
            units2{k} = sprintf( '%02d_%02d', units(k,1), units(k,2));
        end
        units = units2; clear units2
        SUA_chanlabels = unique(units);
        SUA_nchans = length(SUA_chanlabels);
        SUA_nsamps = 1000*diff(sua(1).TimeStamp); % ? ASSUMING here that epoch length will always be cosntant
        SUA_xmin   = sua(1).Event(1,2) - sua(1).Event(2,2); % ASSUMING that stimulus time will not vary across trials relative to t0

        % loop through events and compute spike-density function for MUA and SUA
        SUA_data = zeros(SUA_nchans,SUA_nsamps,EEG.trials);
        MUA_data = zeros(5,SUA_nsamps,EEG.trials);
        for e = 1:EEG.trials

            % get event for this trial
            ent = sua(e).Event;

            % check event code matches the corresponding EEG event code
            code = s.lfp_conds{ent(2,1)==s.lfp_trig};
            if ~strcmp(code,EEG.event(e).type)
                error 'EVENT CODE MISMATCH BETWEEN SUA AND EEG'
            end
            
            % get all spike times for this trial, relative to trial start
            stemp = sua(e).SUA; 
            t0 = ent(1,2); % trial absolute start time (seconds)
            stemp(:,1) = stemp(:,1) - t0; % convert from absolute spike times to times relative to trial start

            % loop through single-units and compute SDF
            for u = 1:SUA_nchans

                % get SDF for this trial/unit pair
                temp = strsplit(SUA_chanlabels{u},'_');
                elecnum = str2double(temp{1});
                unitnum = str2double(temp{2});
                spiketimes = stemp(stemp(:,2)==elecnum & stemp(:,3)==unitnum, 1);

                % compute spike-density function for this unit for this epoch
                dtemp = f_calcSDF(spiketimes, s.sdf_kernelSD, 4); 
                if (length(dtemp) - SUA_nsamps) == 1
                    dtemp = dtemp(1:end-1);
                end
                SUA_data(u,:,e) = dtemp; % store in data matrix
            end % loop through units

            % loop through electrodes and compute MUA SDF
            for c = 1:5
                if any(stemp(:,2)==c)

                    % get spike times for all units for this electrode
                    spiketimes = stemp(stemp(:,2)==c); 

                    % compute spike-density function for this electrode for this epoch
                    dtemp = f_calcSDF(spiketimes, s.sdf_kernelSD, 4); 
                    if (length(dtemp) - SUA_nsamps) == 1
                        dtemp = dtemp(1:end-1);
                    end
                    MUA_data(c,:,e) = dtemp;
                end
            end

        end % loop through trials

        % import SUA to EEGLAB _________________________________________________________________________
        SUA = pop_importdata('dataformat','array', 'data',SUA_data, 'pnts',SUA_nsamps...
            ,'srate',1000,'nbchan',SUA_nchans, 'xmin',SUA_xmin);
        SUA.setname = [ 'merged_SUA ' sesh{sh}];
        
        % copy events
        SUA.event = rmfield(EEG.event,'urevent');
        [~,t0] = min(abs(SUA.times-0));
        t0s = t0 : SUA.pnts : t0 + (SUA.trials-1) * SUA.pnts; 
        for trl = 1:EEG.trials
            SUA.event(trl).latency = t0s(trl);
        end
        
        % define channel (unit) labels 
        clear chanlocs_temp
        for c = 1:SUA_nchans
            chanlocs_temp(c) = struct(...
                'labels', [ 'S' sprintf('%02d',seshnum) '_' SUA_chanlabels{c}],'type','SUA');
        end
        SUA.chanlocs = chanlocs_temp; clear chanlocs_temp
      
        % store original raw spike data in SUA.etc.raw (correcting electrode names/numbers for MK)
        SUA.etc.raw = sua;
        % ________________________________________________________________________________________________

        % import MUA into EEGLAB _________________________________________________________________________
        MUA = pop_importdata('dataformat','array', 'data',MUA_data, 'pnts',SUA_nsamps...
            ,'srate',1000,'nbchan',5, 'xmin',SUA_xmin);
        MUA.setname = [ 'merged_MUA ' sesh{sh}];

        % copy events from SUA
        MUA.event = SUA.event;

        % define channel (unit) labels 
        clear chanlocs_temp
        for c = 1:5
            chanlocs_temp(c) = struct(...
                'labels', ['MUA_' num2str(c)],'type','MUA');
        end
        MUA.chanlocs = chanlocs_temp; clear chanlocs_temp

        % store original raw spike data in MUA.etc.raw (correcting electrode names/numbers for MK)
        MUA.etc.raw = sua;
        % ________________________________________________________________________________________________
    end

    %% baseline-correction 
    labels = {EEG.chanlocs.labels};
    types = {EEG.chanlocs.type};

    % EEG / LFP
    chans2bl = find( strcmp(types,'EEG') | strcmp(types,'LFP') | strcmp(types,'Misc'));
    EEG  = pop_rmbase(EEG, s.blwin_eeg, [], chans2bl ); 

%     % SUA / MUA  %? not doing this for now becuase probably redundant when I do the z-scoring later in VPmonkey_mergeSesh
%                   ? could do z-scoring here also but later means I can re-run it faster
%     if s.sua
%         SUA  = pop_rmbase(SUA, s.blwin_spk, [], 1:SUA.nbchan ); 
%         MUA  = pop_rmbase(MUA, s.blwin_spk, [], 1:MUA.nbchan ); 
%     end

    %% find remaining artifactual pupil epochs and set them all to zero
    if s.eye

        % find epochs which are too flat
        chan = find(strcmp({EEG.chanlocs.labels},'eye_PupilDiameter') );
        pup = squeeze( EEG.data( chan,:,:) );
        pup = pup( 1:16:EEG.pnts ,:); % downsample pupil signal to reduce repeat values from undersampling
        prop_flat = sum( diff(pup)==0 ) / size(pup,1); % compute proportion of timepoints that show change
    %         figure; histogram(prop_flat)
        inds2remove1 = find( prop_flat  > 0.35 );
    
        % find epochs with baselines that are too flat
        puptimes = EEG.times(1:16:EEG.pnts); 
        base = pup(   puptimes >= s.blwin_pup(1) & puptimes <= s.blwin_pup(2)  ,:); clear temp
        prop_flat = sum( diff(base)==0 ) / size(base,1); % compute proportion of timepoints that show change
        inds2remove2 = find( prop_flat  > 0.6 );

        % set bad epochs to zero
        EEG.data(chan,:,inds2remove1) = 0;
        EEG.data(chan,:,inds2remove2) = 0;

% %  ? HERE tried to find trials with too many signal dropouts but results were kinda mixed.. not sure if worth it
%         % find epochs with many larger jumps indicating signal dropouts
%         chan_raw = find(strcmp({EEG.chanlocs.labels},'eye_PupilDiameterRAW') );
%         pup_raw = squeeze( EEG.data( chan_raw,:,:) );
%         % exclude the baseline and bit where visual artifacts are (because these are not removed in the raw signal)
%         puptimes = EEG.times;
%         [~,tstart] = min(abs(puptimes-100)); 
%         [~,tend]   = min(abs(puptimes-3000)); 
%         pup_raw = pup_raw(tstart:tend,:);
%         puptimes = puptimes(tstart:tend);
%         % sliding window mean with very short window to detect fast transients but NOT "flickers"/spikes
%         winlen = round(30 / mean(diff(puptimes))); % time in ms, winlen is samples
%         pup_raw_filt = movmean(pup_raw, winlen); 
%         % find trials with a high percentage of outlying sample-wise difference values
%         diffs = abs(diff(pup_raw_filt));
%         medthresh = median(diffs(diffs>1e-7)); % get Z-threshold of all difference values which are not trivially-zero
%         percout = sum(diffs>(3*medthresh))/size(diffs,1);
%         inds2remove3 = percout>0.15;
%         
%         test1 = pop_select(EEG,'trial',find( inds2remove3)); test1.setname = 'removed';
%         test2 = pop_select(EEG,'trial',find(~inds2remove3)); test2.setname = 'kept';
%         pop_letsplot([test1,test2])
% 
%         % set trial to zeros for easier rejection later
%         EEG.data(chan,:,inds2remove3) = 0;

    end

    %% downsample all channels to 1024  
    dsf = 2;
    EEG = pop_rs_downsample(EEG,dsf);

    % resample SUA to match EEG
    if s.sua
        SUA = pop_resample(SUA,EEG.srate);
        MUA = pop_resample(MUA,EEG.srate);
    end

    %% high-pass filter Pupil channel
    if s.eye
        pupchan = find(strcmp({EEG.chanlocs.labels},'eye_PupilDiameter'));
        PUP = pop_select(EEG,'channel',pupchan);
%         PUP = pop_eegfiltnew( PUP, 'locutoff', 0.05); % maybe a bit too soft
        PUP = pop_eegfiltnew( PUP, 'locutoff', 0.1); % a bit more aggressive to ensure low frequency drift is gone
    %     PUP = pop_rmbase(PUP, s.blwin_pup, [], 1 );   
        %         pop_letsplot([EEG,PUP])
        EEG.data(pupchan,:,:) = PUP.data;
    end

    %% Baseline correction - Pupil channels
    if s.eye
        chans2bl = find( startsWith(labels,'eye_Pupil','IgnoreCase',true) );
        EEG  = pop_rmbase(EEG, s.blwin_pup, [], chans2bl ); 
    end

    %% reject trials that are bad in all datasets and will never be included for any reasons

    % make initial vector of trials to be rejected
    trials2rej = false(size(EEG.data,3),1);

    % exceptions
    if strcmp(subid,'M') && seshnum >= 27 && seshnum <= 29
        trials2rej( strcmp({EEG.event.type},'SOM') ) = true; % reject all somatosensory trials in M for these sessions
%     elseif 
    end
    if seshnum == 7
        trials2rej( strcmp({EEG.event.type},'VIS') ) = true; % reject all visual trials because vis stim stopped working in session S07
    end

%     % convert to fieldtrip and plot data in ft_databrowser 
%     if s.plotData
% 
%         % convert to fieldtrip
%         dataft = eeglab2fieldtrip( EEG, 'raw', 'none' );
% 
%         % plot trials and mark manually which trials to reject
%         trials2rej = VPmonkey_ft_plotTrials(dataft, s.plotYlims2, find(strcmp({EEG.chanlocs.type},'EEG') & ~startsWith({EEG.chanlocs.labels},'EXG')), ...
%             trials2rej, s.chanlay); % eeg only
% %         trials2rej = VPmonkey_ft_plotTrials(dataft, find(strcmp({EEG.chanlocs.type},'EEG')), trials2rej, s.chanlay); % includes EXG chans
% 
%     end
    
    % remove trials   ? only for trials that are bad in all modalities, e.g. the exceptions above
    EEG = pop_select( EEG, 'notrial', find(trials2rej) );
    if s.sua
        SUA = pop_select( SUA, 'notrial', find(trials2rej) );
        MUA = pop_select( MUA, 'notrial', find(trials2rej) );
    end

    %% store original EEG struct for splitting by channel-group in following sections
    EEG_all = EEG;

    %% save EYE tracking data 
    if s.eye

        EYE = pop_select(EEG_all, 'channel', find(strcmp({EEG_all.chanlocs.type},'Eye')) );
        EYE.etc.preproc = s; % save preprocessing settings here
    
        % define save name & path
        savePath = s.savePath;
        saveName = [ 'merged_EYE ' sesh{sh}];
    
        % crop for EEGLAB
        EYE = pop_select(EYE, 'time', s.xlims.lab.EYE);
    
        % save as EEGLAB 
        EYE.setname = saveName;
        pop_saveset(EYE, 'filepath', savePath, 'filename', saveName  );
    
        % downsample for efficiency
        EYE = pop_rs_downsample(EYE,s.dsf);
    
        % crop for lw
        EYE = pop_select(EYE, 'time', s.xlims.lw.EYE);
    
        % replace all NaNs with zeros for easier lw plotting
        EYE.data(isnan(EYE.data)) = 0;
    
        % export to lw
        cd([savePath filesep 'lw'])
        rs_convert_lab2lw_V1( EYE, saveName, [] );
        cd([homedir filesep sesh{sh}]) % go back to session directory
        fprintf('... saved session %s: %s\n', sesh{sh}, [savePath filesep saveName])    
    end

    %% save LFP
    if s.lfp

        LFP = pop_select(EEG_all, 'channel', find(strcmp({EEG_all.chanlocs.type},'LFP')) );
        LFP.etc.preproc = s; % save preprocessing settings here
    
        % define save name & path
        savePath = s.savePath;
        saveName = [ 'merged_LFP ' sesh{sh}];
    
        % crop for EEGLAB
        LFP = pop_select(LFP, 'time', s.xlims.lab.LFP);
    
        % save as EEGLAB 
        LFP.setname = saveName;
        pop_saveset(LFP, 'filepath', savePath, 'filename', saveName  );
    
        % downsample for efficiency
        LFP = pop_rs_downsample(LFP,s.dsf);

        % only take average for lw export & rename to CZ for easier lw plotting with EEG  
        LFP = pop_select(LFP, 'channel', find(strcmp({LFP.chanlocs.labels},'LFP_avg')) );
        LFP.chanlocs(1).labels = 'CZ';
    
        % crop for lw
        LFP = pop_select(LFP, 'time', s.xlims.lw.LFP);
    
        % export to lw
        cd([savePath filesep 'lw'])
        rs_convert_lab2lw_V1( LFP  , saveName, [] );
        cd([homedir filesep sesh{sh}]) % go back to session directory
        fprintf('... saved session %s: %s\n', sesh{sh}, [savePath filesep saveName]) 

    end

    %% save SUA
    if s.sua

        SUA.etc.preproc = s; % save preprocessing settings here
    
        % define save name & path
        savePath = s.savePath;
        saveName = [ 'merged_SUA ' sesh{sh}];
    
%         % crop for EEGLAB   ? this is determined by the range of spikes I get from the epoched spike sorting output
%         SUA = pop_select(SUA, 'time', s.xlims.lab.MUA);
    
        % save as EEGLAB 
        SUA.setname = saveName;
        pop_saveset(SUA, 'filepath', savePath, 'filename', saveName  );
    
        % downsample for efficiency
        SUA = pop_rs_downsample(SUA,s.dsf);
    
        % crop for lw
        SUA = pop_select(SUA, 'time', s.xlims.lw.MUA);
    
        % export to lw
        cd([savePath filesep 'lw'])
        rs_convert_lab2lw_V1( SUA  , saveName, [] );
        cd([homedir filesep sesh{sh}]) % go back to session directory
        fprintf('... saved session %s: %s\n', sesh{sh}, [savePath filesep saveName])

    end

    %% save MUA
    if s.sua

        MUA.etc.preproc = s; % save preprocessing settings here
    
        % define save name & path
        savePath = s.savePath;
        saveName = [ 'merged_MUA ' sesh{sh}];
    
%         % crop for EEGLAB ? this is determined by the range of spikes I get from the epoched spike sorting output
%         MUA = pop_select(MUA, 'time', s.xlims.lab.MUA);
    
        % save as EEGLAB 
        MUA.setname = saveName;
        pop_saveset(MUA, 'filepath', savePath, 'filename', saveName  );

        % replace NaNs with zeros
        MUA.data(isnan(MUA.data)) = 0;
    
        % downsample for efficiency 
        MUA = pop_rs_downsample(MUA,s.dsf);
    
%         % only take average for lw export & rename to CZ for easier lw plotting with EEG  
%         MUA = pop_select(MUA, 'channel', find(strcmp({MUA.chanlocs.labels},'MUA_avg')) );
%         MUA.chanlocs(1).labels = 'CZ';

        % crop for lw
        MUA = pop_select(MUA, 'time', s.xlims.lw.MUA);
    
        % export to lw
        cd([savePath filesep 'lw'])
        rs_convert_lab2lw_V1( MUA  , saveName, [] );
        cd([homedir filesep sesh{sh}]) % go back to session directory
        fprintf('... saved session %s: %s\n', sesh{sh}, [savePath filesep saveName]) 

    end

    %% save MISC channels
    if s.misc

        MISC = pop_select(EEG_all, 'channel', find(strcmp({EEG_all.chanlocs.labels},'photo') | strcmp({EEG_all.chanlocs.labels},'micro')) );
        MISC.etc.preproc = s; % save preprocessing settings here
    
        % define save name & path
        savePath = s.savePath;
        saveName = [ 'merged_MISC ' sesh{sh}];
    
        % crop for EEGLAB
        MISC = pop_select(MISC, 'time', s.xlims.lab.MISC );
    
        % save as EEGLAB 
        MISC.setname = saveName;
        pop_saveset(MISC, 'filepath', savePath, 'filename', saveName  );
    
        % crop for lw
        MISC = pop_select(MISC, 'time', s.xlims.lw.MISC );
    
        % export to lw
        cd([savePath filesep 'lw'])
        rs_convert_lab2lw_V1( MISC , saveName, [] );
        cd([homedir filesep sesh{sh}]) % go back to session directory
        fprintf('... saved session %s: %s\n', sesh{sh}, [savePath filesep saveName])   

    end

    %% save EEG only (pre ICA)
    if s.eeg

        EEG = pop_select(EEG_all, 'channel', find(strcmp({EEG_all.chanlocs.type},'EEG') & ~startsWith({EEG_all.chanlocs.labels},'EXG')) );
        EEG.etc.preproc = s; % save preprocessing settings here
        EEG.etc.numChansInterpolated = numChansInterpolated;
    
        % define save name & path
        savePath = s.savePath;
        if s.asr.useASR, ASRname = num2str(s.asr_thresh); else, ASRname  = 'no'; end
        saveName = [ 'noICA '  'ASR_' ASRname ' reref_' s.refname  ' merged_EEG ' sesh{sh}];
    
        % crop for EEGLAB
        EEG = pop_select(EEG, 'time', s.xlims.lab.EEG);
    
        % save as EEGLAB 
        EEG.setname = saveName;
        pop_saveset(EEG, 'filepath', savePath, 'filename', saveName  );
    
        % downsample for lighter lw storage
        EEG = pop_rs_downsample(EEG,s.dsf);
    
        % crop for lw
        EEG = pop_select(EEG, 'time', s.xlims.lw.EEG );
    
        % export to lw
        cd([savePath filesep 'lw'])
        rs_convert_lab2lw_V1( EEG , saveName, chanlocs_lw );
        cd([homedir filesep sesh{sh}]) % go back to session directory
        fprintf('... saved session %s: %s\n', sesh{sh}, [savePath filesep saveName])     
        fprintf([repmat('_',1,100), '\n\n\n'])

    end

end % end session loop

% cd(s.savePath)
cd([s.savePath filesep 'lw'])
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
   %% OLD linearly interpolate segments of pupil channels with bad value (5) in quality channel
%         % ? could do this also for gaze but the duration may need to be different
        
%         % add new event for quality channel = 5
%         quality = EEG.data( strcmp({EEG.chanlocs.labels},'eye_Quality') ,:);
%         quality = movmean(quality, round(0.5*EEG.srate));
%         [~,lats, W] = findpeaks(quality,'MinPeakDistance', 0.3*EEG.srate);
% %         [~,lats, W] = findpeaks(quality,'MinPeakDistance', 0.3*EEG.srate , 'MinPeakHeight',2.5);
%         widths = W / EEG.srate; clear W
% 
%         % compute window lengths for interpolation using widths of peaks + fixed constant
%         interpwin = widths + 50/1000; 
% 
%         if ~isempty(lats)
%             EEG = eeg_addnewevents(EEG, num2cell(lats), repmat({'badquality'},size(lats)));
%             
%             ctemp = [];
%             ctemp.eventcodes  = {'badquality'};
%             ctemp.timewindow = [ 50 100 ]/1000; 
%             ctemp.channels =  find(startsWith({EEG.chanlocs.labels},'eye_pupil','ignorecase',true));
%             EEG = pop_rs_interpolateShockArtifacts(EEG,ctemp);
%             EEG.event( strcmp({EEG.event.type},'badquality') ) = [];
%         end
%         clear lats quality

        %% OLD filter pupil channels
%         chans_pupil = startsWith({EEG.chanlocs.labels},'eye_pupil','ignorecase',true);
%         EEG_pupil = pop_select( EEG, 'channel',  find(chans_pupil));
%         EEG_pupil = pop_eegfiltnew( EEG_pupil, 'hicutoff', 15); % low-pass
% %         EEG_pupil = pop_eegfiltnew( EEG_pupil, 'hicutoff', 2); % low-pass EXTREME
%         EEG.data(chans_pupil,:) = EEG_pupil.data; % re-insert into data


        %% OLD linearly interpolate segments of pupil channels around peaks at extreme absolute values
%         temp = EEG; % (for superimposing data before and after with letsplot)
% 
%         thresh = -0.06;
% %         mininterpwin = 0.2; % minimum interpolation window, can be expanded for larger widths
%         minpeakdistance = 0.1;
% %         overhang = 50/1000; % at edges of interpolation window, extend this much further and take median of values for interpolation range
% 
%         pup = -EEG.data( strcmp({EEG.chanlocs.labels},'eye_pupilDiameter'),:);
% 
%         [~,lats, W] = findpeaks(pup,'MinPeakDistance', minpeakdistance*EEG.srate , 'MinPeakHeight', thresh ); % ...
% %             ,'MaxPeakWidth', round(maxwidth*EEG.srate) );
%         widths = W / EEG.srate; clear W
% 
%         % compute window lengths for interpolation using widths of peaks + fixed constant
%         interpwin = widths + 50/1000; 
% 
%         if ~isempty(lats)
%             EEG = eeg_addnewevents(EEG, num2cell(lats), repmat({'suddendrop'},size(lats)));
%             
%             ctemp = [];
%             ctemp.eventcodes  = {'suddendrop'};
%             ctemp.timewindow = [interpwin; interpwin]/2; 
%             ctemp.channels = find(startsWith({EEG.chanlocs.labels},'eye_pupil','ignorecase',true));
%             EEG = pop_rs_interpolateShockArtifacts(EEG,ctemp);
%             EEG.event( strcmp({EEG.event.type},'suddendrop') ) = [];
%         end
% 
%         pop_letsplot([temp,EEG])
%         clear temp lats pup

        %% OLD linearly interpolate segments of pupil channels around sudden huge drops in signal
% %         temp = EEG; % (for superimposing data before and after with letsplot)
% 
% % ? this version might get a bit hung up with width measurements because of the lack of signal smoothness
% 
%         thresh = 2;
%         interpwin = 0.3;
%         maxwidth = 0.1;
% 
%         pup = EEG.data( strcmp({EEG.chanlocs.labels},'eye_pupilDiameter'),:);
% %         pup = -(pup - mean(pup,'omitnan')) / std(pup,'omitnan'); %   figure; plot(pup)
%         pup = -(pup - median(pup,'omitnan')) / std(pup,'omitnan'); %   figure; plot(pup)
% 
%         [~,lats] = findpeaks(pup,'MinPeakDistance', interpwin*EEG.srate , 'MinPeakHeight', thresh ...
%             ,'MaxPeakWidth', round(maxwidth*EEG.srate) );
%         if ~isempty(lats)
%             EEG = eeg_addnewevents(EEG, num2cell(lats), repmat({'suddendrop'},size(lats)));
%             
%             ctemp = [];
%             ctemp.eventcodes  = {'suddendrop'};
%             ctemp.timewindow = [interpwin interpwin]/2; 
%             ctemp.channels = find(startsWith({EEG.chanlocs.labels},'eye_pupil','ignorecase',true));
%             EEG = pop_rs_interpolateShockArtifacts(EEG,ctemp);
%             EEG.event( strcmp({EEG.event.type},'suddendrop') ) = [];
%         end
% 
% %         pop_letsplot([temp,EEG])
%         clear temp lats pup

