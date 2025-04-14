% 
% 
% 
%       simple merging, mostly for EEG & pupil, does not account for different LFP/MUA electrodes
% 
%       !! not happy with how files are being handled ... instead get file list from tabl
% 
%%
clc
clearvars
close all

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

saveas = '';
% saveas = '_SPK_thresh_45'

homedir = ['/media/rick/Rick_LabData_3/neuro/iannettilab_ongoing/VPmonkey/data/clean_' subfold saveas];
cd(homedir);

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

% try 
%     eeglab
% catch
%     addlab
% end

%% SETTINGS
s = [];
s.savePath =  homedir; % mkdir(s.savePath)

% s.saveas = 'allSesh'; % name for this set of sessions SHOULD ALWAYS INCLUDE "all" at the start for later string matching
s.saveas = 'allGoodEEG'; % name for this set of sessions SHOULD ALWAYS INCLUDE "all" at the start for later string matching

% epoching 
s.xlims.lw.EEG  = [-0.5 1];  
s.xlims.lw.LFP  = [-0.5 1];  
s.xlims.lw.MUA  = [-0.5 1];  
s.xlims.lw.GAM  = [-0.5 1];  
s.xlims.lw.MISC = [-0.2 0.5];
s.xlims.lw.EYE  = [-0.5 6]; 

s.dsf = 2; % downsample before lw export

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};

% file_types = {'EEG', 'LFP', 'EYE', 'MUA', 'MISC', 'GAM'};
file_types = {'EEG', 'LFP', 'EYE', 'MUA'};
% file_types = {'EYE', 'EEG'};
% file_types = {'GAM'};
% file_types = {'EEG'};

%% LOAD TABLE WITH EEG PREPROC INFO 
tbl_preproc = VPmonkey_load_table_preproc;

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    clear ntrials_all

    %% find session files to merge
    cd(s.savePath)
    files_all = dir; files_all = {files_all.name}; 
    files_all = files_all( endsWith(files_all,[ sub '.set']) & ...
          ~contains(files_all,'all') ...
        & ~contains(files_all,'bin') ...
        & ~contains(files_all,'medsplit') )';

    files_eeg = files_all(  contains(files_all,'merged_EEG') );

    % filter by EEG preproc stage
%     files_eeg = files_eeg(startsWith(files_eeg,'noICA'));
%     files_eeg = files_eeg(startsWith(files_eeg,'cpICA'));
%     files_eeg = files_eeg(startsWith(files_eeg,'yesICA')); 
    files_eeg = files_eeg(startsWith(files_eeg,'LP_30')); 
   
%     % select subset of sessions
%     sel = listdlg('ListString', files_eeg, 'SelectionMode','multiple', 'ListSize',[400,600], 'InitialValue',1:length(files_eeg));
%     files_eeg = files_eeg(sel);

    %% FILTER BY EEG DATA QUALITY

    % filter by subject
    tbl = tbl_preproc;
    tbl = tbl( strcmp(tbl.sub,sub(end)) ,:);
    if length(files_eeg) ~= size(tbl,1), error 'FILE NUMBER MISMATCH!!!'; end

    % filter by bad EEG sessions
    switch s.saveas
        case 'allGoodEEG' % good in at least 2 conditions
            temp = sum( [ tbl.rej_aud tbl.rej_som tbl.rej_vis ], 2 );
            indy = temp <= 1; % i.e. no more than 1 bad condition
            goodsesh = tbl.sesh(indy);

        case 'allPerfEEG' % good in all 3 conditions
            temp = sum( [ tbl.rej_aud tbl.rej_som tbl.rej_vis ], 2 );
            indy = temp == 0; % i.e. no bad conditions
            goodsesh = tbl.sesh(indy);

        case 'allSesh'
            goodsesh = tbl.sesh;
    end

    %% loop through file types
    for ft = 1:length(file_types)
    
        %% get session filenames
        file_type = file_types{ft};
        if strcmp(file_type,'EEG')
            files2load = files_eeg;
        else
            files2load = files_all( endsWith(files_all,[ sub '.set']) & contains(files_all,['merged_' file_type]) )';
        end
        if isempty(files2load)
            error('NO FILES FOUND FOR %s',sub)
        end

        % filter by good sessions
        indy = cellfun( @(x) contains(x,goodsesh), files2load );
        files2load = files2load(indy);
        
        %% load sessions
        clear EEG
        nfiles = length(files2load);
        for f = 1:nfiles
            EEG(f) = pop_loadset('filename',files2load{f},'filepath',s.savePath);
        end
        fprintf('%s ... finished loading\n', file_type)
        
        %% merge sessions
        EEG = pop_mergeset(EEG, 1:nfiles);
        fprintf('%s... finished merging\n', file_type)
    
        %% get number of trials from this dataset for later
        ntrials_all(ft) = EEG.trials;

        %% save merged data in EEGLAB
        if strcmp(file_type,'EEG')
            preproc = strsplit(files_eeg{1});
            preproc = preproc( 1 : (find(strcmp(preproc,'merged_EEG'))-1) );
            preproc = strjoin(preproc,' ');
            filename = [ preproc ' merged_' file_type ' ' s.saveas ' ' sub ];
        else
            filename = [ 'merged_' file_type ' ' s.saveas ' ' sub ];
        end
        filepath = s.savePath;
        pop_saveset(EEG, 'filename',filename,'filepath',filepath);
        fprintf('%s ... finished saving in EEGLAB format\n', file_type)

        %% get t=0 events
        t0 = find(EEG.times==0);    
        t0s = t0 : EEG.pnts : t0 + (EEG.trials-1) * EEG.pnts;  % get t = 0 event markers
        ents = EEG.event( ismember([EEG.event.latency],t0s) );
        ents = struct2table(ents);

        %% set any NaN values to zero for lw export
        EEG.data(isnan(EEG.data)) = 0;
    
        %% split by condition & export to lw
        conds = {'AUD','SOM','VIS'};
        for c = 1:length(conds)
    
            % extract this condition
            cond = conds{c}

            EEG_cond = pop_select(EEG,'trial', find(strcmp(ents.type,cond)) );
            
            % downsample for efficiency
            if ~strcmp( file_types{ft}, 'MISC' )
                EEG_cond = pop_rs_downsample(EEG_cond,s.dsf);
            end

            % if LFP or MUA, just select average for simplicity & rename this to CZ
            if strcmp(file_types{ft},'MUA') || strcmp(file_types{ft},'LFP') || strcmp(file_types{ft},'GAM')
                EEG_cond = pop_select(EEG_cond, 'channel', find(strcmp({EEG_cond.chanlocs.labels},[file_types{ft} '_avg'])) );
                EEG_cond.chanlocs(1).labels = 'CZ';
            end

            % crop for lw
            EEG_cond = pop_select(EEG_cond, 'time', s.xlims.lw.(file_type) );

            % replace all NaNs with zeros for easier lw plotting
            EEG_cond.data(isnan(EEG_cond.data)) = 0;
    
            % export to lw
            cd([filepath filesep 'lw'])
            saveName = [ 'ep ' cond ' ' filename];
            rs_convert_lab2lw_V1( EEG_cond ...
                , saveName, [] );
    
        end
       
        fprintf('%s ... finished exporting each condition to lw\n', file_type)
    
    %% END loop through file types
    end

    if length(unique(ntrials_all)) > 1
        error 'MISMATCHING NUMBER OF TRIALS IN THESE DATASETS'
    end

%% END loop through subs
end

%%
cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

