% 
% 
%       reject epochs from merged data right at the end of processing
% 
%       - Richard Somervail, 2023
%%
clc
clearvars
close all

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

homedir = ['/media/rick/Rick_LabData_3/neuro/iannettilab_ongoing/VPmonkey/data/clean_' subfold];
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

s.saveas = 'allSesh';

% epoching 
s.xlims.lw.EEG  = [-0.5 1];  
s.xlims.lw.LFP  = [-0.5 1];  
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

merged_code = 'allSesh';

subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};

file_types = {'EEG', 'LFP', 'EYE', 'MISC'};
% file_types = {'EEG', 'LFP', 'EYE'};
% file_types = {'LFP', 'EYE', 'MISC'};
% file_types = {'EYE'};
% file_types = {'MUA','LFP'}; % this doesn't really work for LFP / MUA analysis ... better to use other merging script

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    savestr = []; % tracking which modules are used to reject trials

    %% get list of files
    cd(s.savePath)
    files_all = dir; files_all = {files_all.name}; 
    files_all = files_all(  endsWith(files_all,[ sub '.set']) & ...
        contains(files_all,merged_code) & ~contains(files_all,'ar')  )';
    if length(files_all) > length(file_types)
        error 'STILL NEED TO CODE BIT TO HANDLE MULTIPLE COPIES OF THE SAME FILES'
    end
    if isempty(files_all)
        error( 'found no merged-session file for subject %s', sub )
    end

    %% loop through file types
    for ft = 1:length(file_types)
    
        %% get filename for this filetype
        file2load = files_all( endsWith(files_all,[ sub '.set']) &  contains(files_all,['merged_' file_types{ft}]) )';

        data.(file_types{ft}) = pop_loadset('filename',file2load,'filepath',s.savePath);
        
    end

    %% check number of trials in all datasets are matching
    ntrials_all = nan(length(file_types),1);
    for ft = 1:length(file_types)
        ntrials_all(ft) = data.(file_types{ft}).trials;
    end
    if length(unique(ntrials_all)) > 1
        error 'MISMATCHING NUMBER OF TRIALS IN THESE DATASETS'
    end

    %% create master list of epochs for later joint rejection
    temp = data.(file_types{end});
    t0 = find(temp.times==0);    
    t0s = t0 : temp.pnts : t0 + (temp.trials-1) * temp.pnts;  % get t = 0 event markers
    ents = temp.event( ismember([temp.event.latency],t0s) );
    ents = struct2table(ents);
    
    % confirm events correspond uniquely to epochs  
    if size(ents,1) ~= data.(file_types{end}).trials
        error 'number of events doesn''t match number of epochs '; % ? otherwise could select based on latencies (should be monotonic and related to t=0 and trial length)
    end
    ntrials = size(ents,1);

    % add fields for marking a trial as bad
    for ft = 1:length(file_types)
        ents.(['badtrial_' file_types{ft}]) = false(ntrials,1);
    end

    %% AR - auto-reject epochs that have already been marked as bad in pupilometry data
    if isfield(data,'EYE')
        EYE = data.EYE;
        pupchan = strcmp({EYE.chanlocs.labels},'eye_PupilDiameter');
        pup = squeeze(EYE.data( pupchan ,:,:));
    
        % remove any epochs with nan values in pupil channel
        nanpoints = isnan(pup);
        ents.badtrial_EYE(any(nanpoints)) = true;
    
        % remove any epochs with ONLY zero values in pupil channel
        zeropoints = pup == 0;
        ents.badtrial_EYE(sum(zeropoints) == EYE.pnts) = true;
    
        clear EYE
        savestr = [ 'ar_badpupil ' savestr ];
    end

    %% AR - auto-reject epochs based on DEVIANCE from average
    ar_types = {'EYE', 'EEG'}; 
%     ar_types = {'EEG'}; 
    ar_types = {};
    for ft = 1:length(file_types)
        if ismember(file_types{ft}, ar_types)

            % replace zeros with NaNs so that they don't interfere with averaging process
            datatemp = data.(file_types{ft});
            for c = 1:datatemp.nbchan
                chan = squeeze(datatemp.data(c,:,:));
                inds = sum( chan ~= 0 ) == 0;
                datatemp.data(c,:,inds) = nan;
            end

            % config
            ctemp = [];
            switch file_types{ft}
                case 'EYE' % ! re-evaluate this after adding high-pass filter for pupil
                    ctemp.channels = {'eye_PupilDiameter'};
                    ctemp.timewin  = [ 0.1 6 ]; % s 
                    ctemp.thresh = 2;           
                    ctemp.timeprop = 0.2;  
                    ctemp.method = 'median';  
                case 'EEG'
                    ctemp.channels = {data.EEG.chanlocs.labels};
                    ctemp.timewin  = [ 0 0.8 ]; % s 
                    ctemp.thresh = 2;           
                    ctemp.timeprop = 0.2;  
                    ctemp.method = 'median'; 
            end
 
            % run
            [~, badtrials] = pop_rs_autoreject_epochs( datatemp,ctemp );

            % mark bad trials
            ents.(['badtrial_' file_types{ft}])(badtrials) = true;

        end
    end

    savestr = [ 'ar_dev_' ctemp.method '_th' num2str(ctemp.thresh) '_tprop' num2str(100*ctemp.timeprop) ...
           '_' strjoin(ar_types,'_') ' ' savestr];

    %% merge bad trial metrics for each datatype
    badtrials_all = false(ntrials,1);
    for ft = 1:length(file_types)
        badtrials_all = badtrials_all | ents.(['badtrial_' file_types{ft}]);
    end
    trials2keep = find(~badtrials_all);

    %% reject bad trials
    ents = ents(trials2keep,:);
    for ft = 1:length(file_types)
        data.(file_types{ft}) = pop_select(data.(file_types{ft}) ,'trial', trials2keep );
    end

    %% save to EEGLAB
    for ft = 1:length(file_types)
        filename = files_all{ endsWith(files_all,[ sub '.set']) & ...
                              contains(files_all,['merged_' file_types{ft}]) };
        filename = [ savestr filename ];
        pop_saveset( data.(file_types{ft}), 'filename', filename, 'filepath',s.savePath);
        fprintf('saved to:\n%s\n',filename)
    end
    
    %% split by condition & export to lw
    conds = {'AUD','SOM','VIS'};
    for ft = 1:length(file_types)
        EEG = data.(file_types{ft});

        filename = files_all{ endsWith(files_all,[ sub '.set']) & ...
                              contains(files_all,['merged_' file_types{ft}]) };
        filename = [ savestr filename ];

        for c = 1:length(conds)
    
            % extract this condition
            cond = conds{c};
    
            EEG_cond = pop_select(EEG,'trial', find(strcmp(ents.type,cond)) );
            
            % downsample for efficiency
            if ~strcmp( file_types{ft}, 'MISC' )
                EEG_cond = pop_rs_downsample(EEG_cond,s.dsf);
            end
    
            % crop for lw
            EEG_cond = pop_select(EEG_cond, 'time', s.xlims.lw.(file_types{ft}) );
    
            % export to lw
            cd([s.savePath filesep 'lw'])
            saveName = [ 'ep ' cond ' ' filename(1:end-4)];
            rs_convert_lab2lw_V1( EEG_cond ...
                , saveName, [] );
    
        end
   
        fprintf('\n%s ... finished exporting each condition to lw\n\n\n', file_types{ft})
    end

    %% print out details of rejection
    fprintf('\n\n%s:\n',sub)
    fprintf('Rejecting: \t%03d/%03d trials (%.2f%%)\nKeeping: \t%03d/%03d trials (%.2f%%)\n', ...
    sum( badtrials_all), ntrials, 100*sum( badtrials_all)/ntrials, ...
    sum(~badtrials_all), ntrials, 100*sum(~badtrials_all)/ntrials)

    % ! expand this with details about each different type of rejection 
%       (or just by class, e.g. rejecting using EYE or from EEG etc

%% END loop through subs
end

%%
cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

