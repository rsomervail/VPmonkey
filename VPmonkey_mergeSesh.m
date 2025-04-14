% 
% 
% 
%       VPmonkey - merge sessions, either by session or by electrode
%       
%       data = VPmonkey_mergeSesh(cfg)
%           cfg.sub       = 'SubM' or 'SubT'
%           cfg.filetypes = {'EEG','MUA'}; % etc
%           cfg.[filetype]_include = {str}; % for [filetype] (e.g. EEG_include) only include files with string
%           cfg.[filetype]_exclude = {str}; % for [filetype] (e.g. EEG_exclude) only include files without string
% 
%           cfg.byElectrode =       0 - by original session
%                               OR 
%                                   1 - treat each LFP channel as one channel over time & repeat EEG
%                               OR 
%                                   2 - treat each LFP channel as a separate channel, which is nan 
%                                       in trials outside the session it belongs too                            
% 
%           cfg.average   = 'mean' or 'median' or no field to skip averaging
%           cfg.zscore    = cell array of filetypes to be z-scored
%           cfg.zscore_cond = bool, z-score within condition? default = false
%           cfg.zscore_win = time-range for z-scoring (s), [min, max], or struct with cfg.zscore_win.LFP = [min,max] etc
%           cfg.exportlw  = bool, export to lw after merging
%           cfg.filt_elec = e.g. 'Depth > N', where N is number, or can be MED for median value, SPACES are important
%           cfg.filt_preproc  = 'allSesh' or 'allGood' for condition-wise good EEG session selection
%           cfg.filt_eye  = bool, for rejection of sessions with bad single-trial pupil signal (e.g. lots of sensor dropouts)
%           cfg.ar = struct with options for pop_rs_autoreject_epochs          
%           cfg.ar.filetypes = which filetypes to use for pop_rs_autoreject_epochs      
%           cfg.ar_lfp = struct, rejection of trials based on LFP absolute signal and/or std threshold
%               cfg.ar_lfp.thresh_abs = double, absolute amplitude threshold (prior to z-scoring)
%                                       % first column is absolute thresholds, second is max proportion of timepoints above threshold
%               cfg.ar_lfp.maxbad = double, max proportion (0-1) of epochs that can be rejected before electrode is rejected
% 
%           cfg.ar_lfp_SD = reject LFP electrodes by their SD across all epochs
%               cfg.ar_lfp_SD.thresh = threshold in SD over median SD values across epochs
%                                       i.e. if a electrode SD value was 2 SD higher than the median
%                                       (note that the first SD comes from LFP timecourses, 
%                                       while second is the SD of this set of SD values).
%                                       default = 1.5
% 
%           cfg.ar_eye = struct, same as above except with signed threshold (ar_eye.thresh_sign) instead          
%               ? if this works nicely for the LME_EYE_LFP then add to docs and paper properly 
% 
%%
function [data,tbl_depths] = VPmonkey_mergeSesh(cfg)

%% GET PARAMETERS
% homedir = evalin('base','homedir');
homedir = '/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_main';
%
% flds = fields(cfg);
% for k = 1:length(flds)
%     eval([ flds{k} ' = cfg.(flds{k});'  ]);
% end

% defaults
conds = {'AUD','SOM','VIS'};
if ~isfield(cfg, 'zscore')
    cfg.zscore = [];
end
if ~isfield(cfg, 'zscore_cond')
    cfg.zscore_cond = false;
end
if isfield(cfg,'zscore_win')
    if ~isstruct(cfg.zscore_win)
        % if not specified for each datatype listed in cfg.zscore, then repeat for each datatype
        temp = cfg.zscore_win;
        cfg = rmfield(cfg,'zscore_win');
        for k = 1:length(cfg.zscore)
            cfg.zscore_win.(cfg.zscore{k}) = temp;
        end
    end
end
if ~isfield(cfg, 'exportlw')
    cfg.exportlw = false;
end
if ~isfield(cfg, 'mergesesh')
    cfg.mergesesh = true;
end
if ~isfield(cfg,'filt_preproc') % this refers to EEG preproc
    cfg.filt_preproc = 'allSesh';
end
if ~isfield(cfg,'filt_eye')
    cfg.filt_eye = false;
end
if ~isfield(cfg,'filt_elec')
    cfg.filt_elec = '';
end
if ~isfield(cfg,'byElectrode')
    cfg.byElectrode = 0;
end
if isfield(cfg,'ar')
    if ~isfield(cfg.ar,'filetypes')
        cfg.ar.filetypes = cfg.filetypes;
    end
end
if isfield(cfg,'ar_lfp')
    if isempty(cfg.ar_lfp)
        cfg.ar_lfp.thresh_abs = [150,10/100; 300,0]; % good middle-ground (see VPmonkey_LFP_ar_testing for notes)
        cfg.ar_lfp.maxbad     = 0.25; % in practice, could be much higher and still only 1 electrode is removed, from SubT (none for SubM)
    end
end
if ~isfield(cfg,'ar_lfp_SD')
    cfg.ar_lfp_SD.thresh = 1.5;
end

%% GET INFO
ntypes = length(cfg.filetypes);

%% GET PREPROCESSING TABLE
tbl_preproc = VPmonkey_load_table_preproc;
tbl_preproc(  ~strcmp(tbl_preproc.sub, cfg.sub(end)) ,:) = [];
nsesh = length(tbl_preproc.sesh);

%% IF DOING CONDITION-WISE BAD EEG SESSION REMOVAL, THEN REMOVE ANY SESSION WHICH HAS NO GOOD CONDITIONS
if strcmp(cfg.filt_preproc,'allGood')
    temp = sum([ tbl_preproc.rej_aud, tbl_preproc.rej_som, tbl_preproc.rej_vis ],2);
    sesh2remove = temp == 3;
    tbl_preproc(sesh2remove,:) = [];
    nsesh = length(tbl_preproc.sesh);
end

%% PUPIL-WISE BAD SESSION REMOVAL (SINGLE TRIALS)
if cfg.filt_eye
    sesh2remove = tbl_preproc.rej_eye == 1;
    tbl_preproc(sesh2remove,:) = [];
    nsesh = length(tbl_preproc.sesh);
end

%% GET SESSION LIST TABLE
tbl_seshlist = VPmonkey_importSessionlist;

%% GET ELECTRODE DEPTHS TABLE
tbl_depths = VPmonkey_import_electrodeDepths;
tbl_depths(  ~strcmp(tbl_depths.sub, cfg.sub(end))  ,:) = [];

%% FILTER CHOSEN ELECTRODES FROM THE TABLE OF ELECTRODE DEPTHS AND OTHER INFO
tbl_depths(  ~ismember(   tbl_depths.sesh, tbl_preproc.sesh )   ,:) = []; 

% remove NaN electrodes (by default)
tbl_depths = tbl_depths(~isnan(tbl_depths.Depth),:);

% filter any electrode property from electrode table (e.g. Depths, recording site etc)
if ~isempty(cfg.filt_elec)
    % interpret median wild card (MED)
    if contains(cfg.filt_elec,'MED')
        % get chosen field
        field = strsplit(cfg.filt_elec);
        field = field{1};
    
        % check chosen field for any NaN values and remove them
        tbl_depths( isnan(tbl_depths.(field)) ,:) = [];
    
        % compute median
        med = median(tbl_depths.(field));
        cfg.filt_elec = strrep(cfg.filt_elec,'MED',num2str(med));
    end
    % run filter
    eval( sprintf( 'tbl_depths = tbl_depths(tbl_depths.%s,:);', cfg.filt_elec) )
end

%% LOAD FILES
cd(homedir)
files = dir; files = {files.name}';
files = files(endsWith(files,'.set'));

% loop through file-types
for ft = 1:ntypes
    filetype = cfg.filetypes{ft};
    data.(cfg.filetypes{ft}) = [];

    % loop through chosen sessions
    for sh = 1:nsesh
        sesh = tbl_preproc.sesh{sh};

        % filter by file-type, session & subject
        file2load = files( endsWith(files, ['merged_' filetype ' ' sesh '_MT ' cfg.sub '.set']) );      

        % file-type specific filters 
        % include only files containing ALL of the include filters
        if isfield(cfg,[filetype '_include'])
            filt = cfg.([filetype '_include']);
            ind = true(size(file2load));
            for k = 1:length(filt)
                ind = ind & contains( file2load, filt{k} );
            end
            file2load = file2load(ind);
        end
        % exclude any file containing ANY one of the exclude filters
        if isfield(cfg,[filetype '_exclude'])
            filt = cfg.([filetype '_exclude']);
            ind = false(size(file2load));
            for k = 1:length(filt)
                ind = ind | contains( file2load, filt{k} );
            end
            file2load(ind) = []; 
        end

        % check a single file was found
        if length(file2load) > 1
            fprintf('\n\n')
            fprintf('%s\n',file2load{:})
            error('%s: TOO MANY FILES - POSSIBLE FILE AMBIGUITY - FILTERS TOO STRICT OR NOT SPECIFIC ENOUGH?',mfilename)
        elseif isempty(file2load)
            error('%s: NO FILES FOUND',mfilename)
        end

        %% if arranging by electrode method 1 - treating each LFP channel as one channel over time & repeat EEG
        if cfg.byElectrode == 1

                % get any selected electrodes for this session
                tbl = tbl_depths( strcmp(tbl_depths.sesh,sesh) ,:);
                elec2load = tbl.elec;
                
                if ~isempty(elec2load)

                    % IF intracranial data then split by channels 
                    if ~any(strcmp(filetype, {'EEG','EYE','MISC'}))

                        % load electrodes
                        dtemp = pop_loadset('filepath',homedir, 'filename',file2load...
                            ,'loadmode', elec2load);

                        % remove any non-t0 events
                        dtemp = pop_rs_removeNonZeroEvents(dtemp); 

                        % handle selection of good conditions for this session
                        if strcmp(cfg.filt_preproc,'allGood')
                            dtemp = local_remBadConds(dtemp);
                        end

                        % split by electrode
                        dtemp_all = dtemp;
                        for c = 1:length(elec2load)
                            dtemp(c,1) = pop_select(dtemp_all, 'channel',c);
                            dtemp(c,1).chanlocs.labels = filetype;
                        end; clear dtemp_all
        
                    % ELSE then repeat EEG, EYE, MISC data once for each electrode
                    else
                        % load data and repeat
                        dtemp = pop_loadset('filepath',homedir, 'filename',file2load);
                        
                        % remove any non-t0 events
                        dtemp = pop_rs_removeNonZeroEvents(dtemp); 

                        % handle selection of good conditions for this session
                        if strcmp(cfg.filt_preproc,'allGood')
                            dtemp = local_remBadConds(dtemp);
                        end

                        % make a copy for each electrode
                        dtemp = repmat(dtemp,length(elec2load),1);
                    end

                    % store data from this session
                    data.(cfg.filetypes{ft}) = [data.(cfg.filetypes{ft}); dtemp];
                    
                end

        %% ELSE IF arranging by electrode method 2 - treat each LFP channel as a separate channel, which is nan in trials outside the session it belongs too
        elseif cfg.byElectrode == 2

                % get any selected electrodes for this session
                tbl = tbl_depths( strcmp(tbl_depths.sesh,sesh) ,:);
                elec2load = tbl.elec;

                % get any selected electrodes for this session
                tbl = tbl_depths( strcmp(tbl_depths.sesh,sesh) ,:);
                elec2load = tbl.elec;
                
                if ~isempty(elec2load)

                    % IF intracranial data then split by channels 
                    if ~any(strcmp(filetype, {'EEG','EYE','MISC'}))

                        % load electrodes
                        dtemp = pop_loadset('filepath',homedir, 'filename',file2load...
                            ,'loadmode', elec2load);

                        % remove any non-t0 events
                        dtemp = pop_rs_removeNonZeroEvents(dtemp); 

                        % handle selection of good conditions for this session
                        if strcmp(cfg.filt_preproc,'allGood')
                            dtemp = local_remBadConds(dtemp);
                        end

                        % rename electrodes
                        elecnums = find(strcmp(tbl_depths.sesh,sesh)); % get electrode numbers from table of electrodes
                        for c = 1:dtemp.nbchan
                            dtemp.chanlocs(c).labels = sprintf('LFP_%03d',elecnums(c));
                        end
            
                    % ELSE then repeat EEG, EYE, MISC data once for each electrode
                    else
                        % load data and repeat
                        dtemp = pop_loadset('filepath',homedir, 'filename',file2load);
                        
                        % remove any non-t0 events
                        dtemp = pop_rs_removeNonZeroEvents(dtemp); 

                        % handle selection of good conditions for this session
                        if strcmp(cfg.filt_preproc,'allGood')
                            dtemp = local_remBadConds(dtemp);
                        end

                    end

                    % store data from this session
                    data.(cfg.filetypes{ft}) = [data.(cfg.filetypes{ft}); dtemp];
                    
                end % if empty elec2load

        %% ELSE IF NOT PROCESSING BY ELECTRODE
        else
            % load data for this session
            dtemp = pop_loadset('filepath',homedir, 'filename',file2load);

            % remove any non-t0 events
            dtemp = pop_rs_removeNonZeroEvents(dtemp);

            % handle selection of good conditions for this session
            if strcmp(cfg.filt_preproc,'allGood')
                dtemp = local_remBadConds(dtemp);
            end

            % store data from this session
            if isempty(data.(cfg.filetypes{ft}))
                data.(cfg.filetypes{ft}) = dtemp;
            else
                data.(cfg.filetypes{ft})(sh) = dtemp;
            end

        end % end byElectrode check

    end % sesh loop

    % get number of files in case byElectrode or byUnit was used
    nfiles = length(data.(cfg.filetypes{1}));
  
    %% IF EYE, then remove all channels except pupil diameter
    if strcmp(cfg.filetypes{ft},'EYE')
        for f = 1:nfiles
        
            data.(cfg.filetypes{ft})(f) = pop_select(data.(cfg.filetypes{ft})(f), ...
                'channel', 8);

            if ~strcmp({data.(cfg.filetypes{ft})(f).chanlocs.labels},'eye_PupilDiameter')
                error 'GOT THE WRONG EYE CHANNEL ...'
            end
%             data.(cfg.filetypes{ft})(f) = pop_select(data.(cfg.filetypes{ft})(f), ...
%                 'channel', find(strcmp({data.(cfg.filetypes{ft})(f).chanlocs.labels},'eye_PupilDiameter') ));

        end % session loop
    end

end % filetype loop

%% IF USING BYELECTRODE LOADING AND LFP WAS REQUESTED, REJECT DEVIANT ELECTRODES
if cfg.byElectrode == 1 && any(strcmp(cfg.filetypes,'LFP')) 

    % loop through electrodes
    sds = nan(nsesh,1);
    for sh = 1:nfiles
        % compute SD for the following outlier detection
        dtemp = data.LFP(sh).data;
        sds(sh) = std(dtemp,1,2:3);
    end

    % z-score the SD values for each electrode
%     sds = (sds-mean(sds))/std(sds,1);
    sds = (sds-median(sds))/std(sds,1);

    % remove electrodes exceeding a certain threshold
    ind = find(sds > cfg.ar_lfp_SD.thresh);
    fprintf( 'VPmonkey_mergeSesh: %d/%d electrodes rejected due to high SD (%d%%)\n', ...
        length(ind), length(sds), round(100*length(ind)/length(sds))  )  
    if ~isempty(ind)
        for k = 1:length(ind)
            str = sprintf('%s-%d', tbl_depths.sesh{ind(k)}, tbl_depths.elec(ind(k)));
            fprintf('rejecting electrode: %s\n',str)
        end
    end

    % loop through data and remove bad electrodes
    if ~isempty(ind)
        for k = 1:length(ind)
            for ft = 1:length(cfg.filetypes)
                filetype = cfg.filetypes{ft};
                data.(filetype)(ind(k)).data = [];
                data.(filetype)(ind(k)).chanlocs = [];
                data.(filetype)(ind(k)).nbchan = 0;
            end
        end
    end
    % ? corresponding tbl_depths electrodes are removed later because the dataset will be empty

end

%% IF NOT USING BYELECTRODE LOADING AND LFP WAS REQUESTED, MERGE ELECTRODES AND REJECT DEVIANTS
% ? includes built-in rejection of electrodes with deviant SDs
if cfg.byElectrode == 0 && any(strcmp(cfg.filetypes,'LFP'))

    % loop through sessions
    sds = cell(nsesh,1);
    for sh = 1:nsesh
        dtemp = data.LFP(sh).data;
        % get the electrodes which have non-nan depth values
        ind = strfind(data.LFP(sh).setname,'_MT Sub');
        sesh = data.LFP(sh).setname(ind-3:ind-1);
        ind = tbl_depths{strcmp(tbl_depths.sesh,sesh),"elec"};
        dtemp = dtemp(ind,:,:);
        % compute SD of each remaining channel for the following outlier detection
        sds{sh} = std(dtemp,1,2:3);
    end

    % concatenate and z-score the SD values for each electrode
    sds = cat(1,sds{:});
%     sds = (sds-mean(sds))/std(sds,1);
    sds = (sds-median(sds))/std(sds,1);

    % remove electrodes exceeding a certain threshold
    ind = find(sds > cfg.ar_lfp_SD.thresh);
    fprintf( 'VPmonkey_mergeSesh: %d/%d electrodes rejected due to high SD (%d%%)\n', ...
        length(ind), length(sds), round(100*length(ind)/length(sds))  ) 
    if ~isempty(ind)
        for k = 1:length(ind)
            str = sprintf('%s-%d', tbl_depths.sesh{ind(k)}, tbl_depths.elec(ind(k)));
            fprintf('rejecting electrode: %s\n',str)
        end
    end
    tbl_depths(ind,:) = [];

    % loop through sessions and remove those electrodes, then replace all data with average of remaining channels
    for sh = 1:nsesh
        dtemp = data.LFP(sh).data;
        % get the electrodes which have non-nan depth values and **have not been removed as deviant in the previous step**
        indy = strfind(data.LFP(sh).setname,'_MT Sub');
        sesh = data.LFP(sh).setname(indy-3:indy-1);
        elecs2keep = tbl_depths{strcmp(tbl_depths.sesh,sesh),"elec"};
        % if any electrodes remain, then take average and reinsert into data
        if ~isempty(elecs2keep)
            dtemp = dtemp(elecs2keep,:,:);
            dtemp = mean(dtemp,1);
            data.LFP(sh).data = dtemp;
            data.LFP(sh).chanlocs = data.LFP(sh).chanlocs(strcmp({data.LFP(sh).chanlocs.labels},'LFP_avg'));
            data.LFP(sh).nbchan = 1;
        else
            % if no electrodes remain for this session, then remove all data from other file-types too
            for ft = 1:length(cfg.filetypes)
                filetype = cfg.filetypes{ft};
                data.(filetype)(sh).data = [];
                data.(filetype)(sh).chanlocs = [];
                data.(filetype)(sh).nbchan = 0;
            end
        end
    end
end

%% TRIAL REJECTION STEP 1 - DEFINE VECTOR OF ALL TRIALS TO ALLOW COMMON REMOVAL OF BAD TRIALS
trials2keep1 = cell(nfiles,1);
for f = 1:nfiles 
    trials2keep1{f} = true(data.(cfg.filetypes{1})(f).trials,1);
end

%% TRIAL REJECTION STEP 1 - IF EYE DATA ARE LOADED, REMOVE ALL TRIALS CONTAINING NANS OR CONTAINING ONLY ZEROES
if any(contains(cfg.filetypes,'EYE'))
    for f = 1:nfiles
        dtemp = data.EYE(f);
        nantrials  = squeeze(  any(isnan(dtemp.data)) );
        zerotrials = squeeze( ~any(dtemp.data~=0)     ); 
        trials2keep1{f}(nantrials | zerotrials) = false;
    end


    % print number of rejected epochs (up to and including this step) 
    trials2keep1_EYE = cat(1,trials2keep1{:});
    fprintf( 'VPmonkey_mergeSesh: STEP 1 trial rejection (EYE DATA)\n')
    fprintf( 'VPmonkey_mergeSesh: %d/%d epochs rejected (up to and including this step)  (%d%%)\n', ...
        sum(~trials2keep1_EYE), length(trials2keep1_EYE), round(100*sum(~trials2keep1_EYE)/length(trials2keep1_EYE))  )  
    fprintf( 'VPmonkey_mergeSesh: %d epochs remaining (%d%%)\n', ...
        sum(trials2keep1_EYE), round(100*sum(trials2keep1_EYE)/length(trials2keep1_EYE)) )  

end

%% TRIAL REJECTION STEP 1 - IF EYE DATA ARE LOADED & THIS SETTING IS ENABLED, APPLY **NEGATIVE** SIGNED THRESHOLD
if any(contains(cfg.filetypes,'EYE')) && isfield(cfg,'ar_eye')
    thresh_sign = cfg.ar_eye.thresh_sign;
    for f = 1:nfiles
        dtemp = pop_select(data.EYE(f),'time',cfg.ar_eye.timewin);
        ttemp = dtemp.times/1000;
        dtemp = squeeze(dtemp.data);
        
        % apply each threshold that is used
        badtrials_eye = false(1,size(dtemp,2));
        for k = 1:size(thresh_sign,1)
            badtrials_eye = badtrials_eye | ...
            sum(  dtemp < thresh_sign(k,1)) / size(dtemp,1) ... % 
             > thresh_sign(k,2); % percentage of timepoints threshold
        end

        % if more than N% of trials are removed then remove whole electrode
        if sum(badtrials_eye)/length(badtrials_eye) > cfg.ar_eye.maxbad
            badtrials_eye(:) = true;
        end

% %         % SANITY CHECK PLOT
% %         if any(badtrials_lfp)
% %         if f > 20
%             figure(332);
%             clf
%             for k = 1:size(dtemp,2)
%                 if badtrials_eye(k)
%                     subplot(1,2,1); plot(ttemp, dtemp(:,k),'Color','r','LineWidth',1); hold on
%                 else
%                     subplot(1,2,2); plot(ttemp, dtemp(:,k),'Color','k','LineWidth',1); hold on
%                 end
%             end
%             subplot(1,2,1); ylim([-8,8]/100); title([cfg.sub ' ' num2str(f) '/' num2str(nfiles)]); 
%             subplot(1,2,2); ylim([-8,8]/100); 
% 
%             true; % breakpoint
% %         end

        %
        trials2keep1{f}(badtrials_eye) = false;
    end

end

%% TRIAL REJECTION STEP 1 - IF LFP DATA ARE LOADED, REMOVE ALL TRIALS WITH SIGNAL GREATER THAN THRESHOLD
if any(contains(cfg.filetypes,'LFP')) && isfield(cfg,'ar_lfp')
    thresh_abs = cfg.ar_lfp.thresh_abs;
    for f = 1:nfiles
        dtemp = squeeze(data.LFP(f).data);

        % apply each threshold that is used
        badtrials_lfp = false(1,size(dtemp,2));
        for k = 1:size(thresh_abs,1)
            badtrials_lfp = badtrials_lfp | ...
            sum(  dtemp < -thresh_abs(k,1) | dtemp > thresh_abs(k,1)  ) / size(dtemp,1) ... % 
             > thresh_abs(k,2); % percentage of timepoints threshold
        end

        % if more than N% of trials are removed then remove whole electrode
        if sum(badtrials_lfp)/length(badtrials_lfp) > cfg.ar_lfp.maxbad
            badtrials_lfp(:) = true;
        end

% %         % SANITY CHECK PLOT
% %         if any(badtrials_lfp)
% %         if f > 20
%             figure(332);
%             clf
%             for k = 1:size(dtemp,2)
%                 if badtrials_lfp(k)
%                     subplot(1,2,1); plot(data.LFP(1).times/1000, dtemp(:,k),'Color','r','LineWidth',1); hold on
%                 else
%                     subplot(1,2,2); plot(data.LFP(1).times/1000, dtemp(:,k),'Color','k','LineWidth',1); hold on
%                 end
%             end
%             subplot(1,2,1); ylim([-600,600]); title([cfg.sub ' ' num2str(f) '/' num2str(nfiles)]); 
%             subplot(1,2,2); ylim([-600,600]); 
% 
%             true; % breakpoint
% %         end

        %
        trials2keep1{f}(badtrials_lfp) = false;
    end

    % print number of rejected epochs (up to and including this step) 
    trials2keep1_LFP = cat(1,trials2keep1{:});
    fprintf( 'VPmonkey_mergeSesh: STEP 1 trial rejection (LFP DATA)\n')
    fprintf( 'VPmonkey_mergeSesh: %d/%d epochs rejected (up to and including this step) (%d%%)\n', ...
        sum(~trials2keep1_LFP), length(trials2keep1_LFP), round(100*sum(~trials2keep1_LFP)/length(trials2keep1_LFP))  )  
    fprintf( 'VPmonkey_mergeSesh: %d epochs remaining (%d%%)\n', ...
        sum(trials2keep1_LFP), round(100*sum(trials2keep1_LFP)/length(trials2keep1_LFP)) )  

end

%% TRIAL REJECTION STEP 1 - REJECT ALL TRIALS THAT HAVE BEEN REMOVED SO FAR (BEFORE AUTOREJECTION)
for ft = 1:ntypes
    filetype = cfg.filetypes{ft};
    for f = 1:nfiles
        if any(~trials2keep1{f}) 
            if all(~trials2keep1{f})
                data.(filetype)(f).data = []; % remove data manually so pop_select doesn't complain
            else
                data.(filetype)(f) = pop_select(data.(filetype)(f), 'trial', find(trials2keep1{f}));
            end
        end
    end
end

%% TRIAL REJECTION STEP 2 - RE-DEFINE VECTOR OF ALL REMAINING TRIALS TO ALLOW COMMON REMOVAL OF BAD TRIALS
trials2keep2 = cell(nfiles,1);
for f = 1:nfiles 
    trials2keep2{f} = true(data.(cfg.filetypes{1})(f).trials,1);
end

%% TRIAL REJECTION STEP 2 - FIND BAD TRIALS USING pop_rs_autoreject_epochs
if isfield(cfg,'ar')

    % loop through file-types
    for ft = 1:length(cfg.ar.filetypes)
        filetype = cfg.ar.filetypes{ft};

        % loop through chosen sessions
        for f = 1:nfiles
            if ~isempty(data.(filetype)(f).data) % if no trials for this file remain from previous step then skip this step 
                [~,badtrials] = pop_rs_autoreject_epochs(data.(filetype)(f), cfg.ar);
                trials2keep2{f}(badtrials) = false;
            end
        end
    end
end

%% TRIAL REJECTION STEP 2 - REJECT THE TRIALS IDENTIFIED BY pop_rs_autoreject_epochs
for ft = 1:ntypes
    filetype = cfg.filetypes{ft};
    for f = 1:nfiles
        if any(~trials2keep2{f})
            data.(filetype)(f) = pop_select(data.(filetype)(f), 'trial', find(trials2keep2{f}));
        end
    end
end

%% REMOVE ALL EMPTY SESSIONS/ELECTRODES
% ? assuming empty sessions are the same in all datasets by now
files2remove = arrayfun(@(x) isempty(x.data), data.(cfg.filetypes{1}));
if any(files2remove)
    nfiles = sum(~files2remove);
    if cfg.byElectrode == 0
        nsesh = nfiles;
    end

    % remove empty datasets
    for ft = 1:ntypes
        data.(cfg.filetypes{ft})(files2remove) = [];
    end
    
    % handle electrode depth table
    if cfg.byElectrode == 1  
        tbl_depths(files2remove,:) = [];
        % ? could edit seshlist here too
    elseif cfg.byElectrode == 0

        % removing here any remaining tbl_depths electrodes that were not already removed in the LFP SD step
        sesh2remove = arrayfun(@(x) char(x), tbl_seshlist.session, 'UniformOutput',false);
        sesh2remove = sesh2remove(files2remove);
        elecs2remove = ismember(tbl_depths.sesh, sesh2remove);
        if any(elecs2remove)
            tbl_depths(elecs2remove,:) = [];
        end

    elseif cfg.byElectrode == 2
        error 'should probably just scrap byElectrode == 2'
    end
end

%% SANITY CHECK THAT NUMBER OF TRIALS AND EVENTS MATCH BETWEEN FILE-TYPES

% check trial number is the same across filetypes
ntrials_final = nan(1,ntypes);
for ft = 1:ntypes
    ntrials_final(ft) = data.(cfg.filetypes{ft}).trials;
end
if length(unique(ntrials_final)) ~= 1
    error 'TRIAL NUMBER MISMATCH BETWEEN DATASETS'
else
    ntrials_final = ntrials_final(1);
end

% check that event codes match across filetypes
for f = 1:nfiles
    for ft1 = 1:ntypes
        for ft2 = 1:ntypes
            entcodes1 = {data.(cfg.filetypes{ft1})(f).event.type};
            entcodes2 = {data.(cfg.filetypes{ft2})(f).event.type};
            if ~isequal( entcodes1,entcodes2 )
                error 'EVENT CODE MISMATCH'
            end
        end
    end
end


%% PRINT OUT RESULTS OF TRIAL REJECTION
% trials2keep1
trials2keep1 = cat(1,trials2keep1{:});
fprintf( 'VPmonkey_mergeSesh: STEP 1 trial rejection (file-type specific rejection)\n')
fprintf( 'VPmonkey_mergeSesh: %d/%d epochs rejected (%d%%)\n', ...
    sum(~trials2keep1), length(trials2keep1), round(100*sum(~trials2keep1)/length(trials2keep1))  )  
fprintf( 'VPmonkey_mergeSesh: %d epochs remaining (%d%%)\n', ...
    sum(trials2keep1), round(100*sum(trials2keep1)/length(trials2keep1)) )  

% trials2keep2
trials2keep2 = cat(1,trials2keep2{:});
fprintf( 'VPmonkey_mergeSesh: STEP 2 trial rejection (automatic trial rejection by deviance from avg)\n')
fprintf( 'VPmonkey_mergeSesh: %d/%d epochs rejected (%d%%)\n', ...
    sum(~trials2keep2), length(trials2keep2), round(100*sum(~trials2keep2)/length(trials2keep2))  )  
fprintf( 'VPmonkey_mergeSesh: %d epochs remaining (%d%%)\n', ...
    sum(trials2keep2), round(100*sum(trials2keep2)/length(trials2keep2)) )  


%% Z-SCORE DATA WITHIN SESSION

% loop through file-types
for ft = 1:ntypes
    filetype = cfg.filetypes{ft};

    if ismember(filetype,cfg.zscore)

        % loop through chosen sessions (or electrodes or units)
        for f = 1:nfiles

            % get this session
            dtemp = double(data.(filetype)(f).data);
    
            % Z-SCORE WITHIN SESSION
            if ~cfg.zscore_cond
       
                % handle z-score window
                if isfield(cfg,'zscore_win')
                    times = data.(filetype)(f).times/1000;
                    [~,indy(1)] = min(abs(times - cfg.zscore_win.(filetype)(1))); 
                    [~,indy(2)] = min(abs(times - cfg.zscore_win.(filetype)(2))); 
                    zwin = dtemp(:,indy(1):indy(2),:);
                else
                    zwin = dtemp;
                end
    
                % z-score
                if strcmp(filetype,'EEG')
                    sd = std( zwin ,[],'all'); 
                    if sd ~= 0 % to prevent division by zero
                        dtemp = dtemp / sd; 
                    end
                else
                    sd = std( zwin ,[],[2,3]);  
                    for c = 1:size(dtemp,1)
                        if sd(c) ~= 0 % to prevent division by zero
                            dtemp(c,:,:) = dtemp(c,:,:) / sd(c); 
                        end
                    end
                end

            % Z-SCORE WITHIN CONDITION WITHIN SESSION
            else

                % loop through conds
                for cond = 1:length(conds)

                    % get the trials for this condition
                    condinds   = strcmp({data.(filetype)(f).event.type},conds{cond});
                    dtemp_cond = dtemp(:,:,condinds);
           
                    % handle z-score window
                    if isfield(cfg,'zscore_win')
                        times = data.(filetype)(f).times/1000;
                        [~,indy(1)] = min(abs(times - cfg.zscore_win.(filetype)(1))); 
                        [~,indy(2)] = min(abs(times - cfg.zscore_win.(filetype)(2))); 
                        zwin = dtemp_cond(:,indy(1):indy(2),:);
                    else
                        zwin = dtemp_cond;
                    end
        
                    % z-score
                    if strcmp(filetype,'EEG')
                        sd = std( zwin ,[],'all'); 
                        if sd ~= 0 % to prevent division by zero
                            dtemp_cond = dtemp_cond / sd; 
                        end
                    else
                        sd = std( zwin ,[],[2,3]);  
                        for c = 1:size(dtemp_cond,1)
                            if sd(c) ~= 0 % to prevent division by zero
                                dtemp_cond(c,:,:) = dtemp_cond(c,:,:) / sd(c); 
                            end
                        end
                    end

                    % insert the trials from this condition back into the data
                    dtemp(:,:,condinds) = dtemp_cond;
        
                end % condition loop

            end % if z-score within cond

            data.(filetype)(f).data = single(dtemp);
    
        end % sesh loop

        fprintf('%s: ... finished within-session z-scoring for %s\n',mfilename, filetype)
    end

end % filetype loop

%% AVERAGE TRIALS WITHIN SESSION
if cfg.average

    % loop through file-types
    for ft = 1:ntypes
        filetype = cfg.filetypes{ft};
    
        % WITHIN-SESSION AVERAGE
        data.(filetype) = pop_rs_average(data.(filetype), true);

        fprintf('%s: ... finished averaging within-session for %s\n',mfilename, filetype)
    end % filetype loop

    % handle electrode depths table (matching it to averages becomes unclear 
    %                                 because of different conditions per sessions)
    % loop through sessions within any filetype
    if cfg.byElectrode==1 % ? maybe would be useful too with version 2 if I ever use it?
        inds = cell(nfiles,1);
        for k = 1:nfiles
            inds{k} = repmat(k,data.(filetype)(k).trials,1);
        end
        inds = cat(1,inds{:});
        tbl_depths = tbl_depths(inds,:); clear inds
    end

end % 

%% MERGE ACROSS SESSIONS
if cfg.mergesesh

    %% IF USING BY ELECTRODE METHOD 2 THEN ADD MISSING ELECTRODES CHANNELS AND FILLING MISSING TRIALS WITH NANS)

    % loop through electrode types
    if cfg.byElectrode == 2
        electypes = find(ismember(cfg.filetypes, {'LFP','MUA'}));
        for ft = 1:length(electypes)
            filetype = cfg.filetypes{electypes(ft)};

            % make chanlocs
            templocs = struct;
            for c = 1:size(tbl_depths,1)
                templocs(c).labels = sprintf('%s_%03d',filetype,c);
                templocs(c).type = filetype;
            end

            % loop through files
            for f = 1:nfiles
                % make nan matrix with the size of all channels, missing or present
                dold = data.(filetype)(f).data;
                dnew = nan(size(tbl_depths,1),size(dold,2),size(dold,3));

                % insert present channels into full matrix
                chans_present = ismember({templocs.labels},{data.(filetype)(f).chanlocs.labels});
                dnew(chans_present,:,:) = dold;
                data.(filetype)(f).data = dnew;

                % update dataset info
                data.(filetype)(f).chanlocs = templocs;
                data.(filetype)(f).nbchan   = length(templocs);
            end
        end
    end
        
    %% loop through file-types and MERGE
    for ft = 1:ntypes
        filetype = cfg.filetypes{ft};

            data.(filetype) = pop_rs_mergeEpochs(data.(filetype));
%             data.(filetype) = pop_mergeset(data.(filetype), 1:length(data.(filetype)));

            % remove ICA info which no longer applies to the merged data
            data.(filetype).icaact = [];
            data.(filetype).icawinv = [];
            data.(filetype).icaweights = [];
            data.(filetype).icasphere = [];
            %? chansinds as well
    end

    %% RENAME FILE TO REFLECT MERGED DATA
    % define output dataset name
    if isfield(cfg,[filetype '_include'])
        include = strjoin(cfg.([filetype '_include']));
    else 
        include = '';
    end
    if isfield(cfg,[filetype '_exclude'])
        exclude = strjoin(cfg.([filetype '_exclude']));
    else 
        exclude = '';
    end
    if cfg.byElectrode == 1
        bystr = 'byElec_1';
    elseif cfg.byElectrode == 2
        bystr = 'byElec_2';
    else
        bystr = 'bySesh';
    end 
    if cfg.average
        avgstr = 'avg';
    else
        avgstr = 'trl';
    end 
    if ~isempty(cfg.zscore)
        zstr = 'zYES';
    else
        zstr = 'zNO';
    end 
    if isfield(cfg,'autoar')
        ar = cfg.ar;
        arstr = ['ar_', num2str(ar.thresh)];
    %                 arstr = strjoin({ 'ar_', num2str(ar.thresh)) ...
    % %                    ['t_' num2str(ar.timeprop)]  },'_');
    else
        arstr = 'ar_NO';
    end
    saveName = strjoin( { avgstr arstr zstr ...
      include, exclude, ...
      bystr,  ...
      'merged_$FILETYPE$', cfg.filt_preproc 'group' cfg.sub });
    saveName = strrep(saveName,'   ', ' '); % trim excess whitespaces
    saveName = strrep(saveName,'  ', ' ');
    
    % loop through file-types and rename merged data
    for ft = 1:ntypes
        filetype = cfg.filetypes{ft};
        data.(filetype).setname = strrep(saveName,'$FILETYPE$',filetype);
    end

    %% EXPORT TO LW
    if cfg.exportlw
    
        cd([homedir '/lw/'])
    
        % load good channel locations
        load([getRoot 'EEG Configurations' filesep 'monkeyEEG_LW_chanlocs.mat'],'chanlocs_lw');
    
        % loop through file-types
        for ft = 1:ntypes
            filetype = cfg.filetypes{ft};
    
            EEG = data.(filetype);

            % rename for filetype
            saveName_type = strrep(saveName,'$FILETYPE$',filetype);
    
            % rename single channel file-types to Cz for letswave plotting
            if ~strcmp(filetype,'EEG')
                EEG.chanlocs(1).labels = 'Cz';
            end
    
            % loop through conditions and export
            conds = unique({EEG.event.type});
            for c = 1:length(conds)
        
                % extract this condition
                cond = conds{c};
                EEG_cond = pop_select(EEG,'trial', find(strcmp({EEG.event.type},cond)) );
                
        %         % replace all NaNs with zeros for easier lw plotting
        %         EEG_cond.data(isnan(EEG_cond.data)) = 0;
          
                % export to lw
                if strcmp(filetype,'EEG')
                    rs_convert_lab2lw_V1( EEG_cond ...
                        , ['ep ' cond ' ' saveName_type], chanlocs_lw );
                else
                    rs_convert_lab2lw_V1( EEG_cond ...
                        , ['ep ' cond ' ' saveName_type], [] );
                end
     
            end
        end 
        cd(homedir)
    end

%% END IF MERGE SESH
end

%% SUBFUNCTIONS
function dtemp = local_remBadConds(dtemp)

    % find which conditions that exist in this dataset should be stored
    seshconds = unique({dtemp.event.type});
%     if length(seshconds) ~= 3
%         error '!! must account for missing condition existing in this session!' 
%     end

    % out of the conditions inside this dataset, which conditions should be kept and which rejected?
    indlocal = false(1,length(seshconds));
    for clocal = 1:length(seshconds)
        indlocal(clocal) = ~tbl_preproc.(['rej_' lower(seshconds{clocal})])(sh);
    end
    conds2take = seshconds(indlocal);
    
    % get conditions for each trial
    trialconds = {dtemp.event.type};
    if length(trialconds) ~= dtemp.trials
        error 'too many trials ... non-t0 events found? ... these should be removed from data or handled here'
    end

    % select only trials with the correct conditions
    trials2take = find(ismember(trialconds,conds2take)); 
    if isempty(trials2take)
        error 'NO TRIALS FOUND WHILE DOING CONDITION-WISE BAD EEG REJECTION - WHY?'
    else
        if ~isequal(trials2take,1:dtemp.trials) % no need to select trials if all trials would be selected anyway
            dtemp = pop_select(dtemp,'trial', trials2take);
        end
    end

end


%% FUNCTION END
end

