% 
% 
% 
%       file merging that bins by electrode depth 
% 
%       !! I THINK BEST APPROACH TO BINNING IS STILL TO FIND EDGES THAT MAXIMISE CLUSTERING
%           METRIC, LIKE EVEN IF THINGS AREN'T THAT CLEAR AT LEAST IT IS THE BEST IT CAN BE ...
%           MOST TARGETTED TO WHAT I'M ASKING ...
% 
% 
%         !! SHOULD CHANGE FILENAMES TO BIN_1_3 BIN_3_3 ETC
% 
%           !! got more trials than expected per bin & per condition - does this make sense? maybe cause epochs are 77chan x trial pairs?
%             
%           ? possibly outdated notes here:
%           ? right now I'm taking the averages rather than concatenating trials regardless of electrode
%             could take trials instead, or as well ... and do trial by trial correlations and such
%           ! probs nicer to do this but maybe need some way to keep track of different electrodes (e.g. EEG.event)
%               so that I can later take averages within each and do stats across these averages
% 
% 
%%
clc
clearvars
close all

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

saveas = '';
% saveas = '_SPK_thresh_45';

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

s.saveas = 'allSesh'; % name for this set of sessions SHOULD ALWAYS INCLUDE "all" at the start for later string matching

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

s.conds = {'AUD','SOM', 'VIS'};

subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};

% file_types = {'GAM'};
% file_types = {'MUA'};
file_types = {'LFP'};
% file_types = {'MUA','LFP','EEG'}; ?? to extract EEG also I'll need to add a bit that will take the corresponding EEG for each LFP/MUA within each bin
% file_types = {'MUA','LFP','EEG','EYE'};
% file_types = {'MUA','LFP','EEG','EYE','MISC'};

% depth metric to use for binning
binmetric = 'Depth';
% binmetric = 'DepthRelDura';
% binmetric = 'DepthRelTOA';

% method to use for binning 
binmethod = 'abs';  % bin by absolute values, results in different numbers within each bin
% binmethod = 'rank'; % bin by rank, results in equal numbers within each bin i.e. with nbins = 2 this is median split
                        % ? for equal statistical power between bins

% number of bins to use
nbins = 3;

%% import intracranial information table
tbl_depths = VPmonkey_import_electrodeDepths;

%% select sessions (assumes the necessary files exist)
sesh = tbl_depths.sesh;
% sel = listdlg('ListString',sesh,'SelectionMode','multiple');
% sesh = sesh(sel);
nfiles = length(sesh);

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    clear ntrials_all

    %% define depth-bins given the selected sessions & subject
    % get all depths for the chosen bin metric
    inds = cellfun( @(x) find(strcmp(x, tbl_depths.sesh))  , sesh );
    deps_all = [];
    for k = 1:5
        temp = tbl_depths.([binmetric '_' sub(end) num2str(k)]);
        deps_all = [ deps_all; temp(inds) ]; %#ok<*AGROW> 
    end
    deps_all(isnan(deps_all)) = []; %#ok<*SAGROW> 
    deps_all = sort(deps_all); % ? this vector should not be used for indexing later channels anyway
                                   %   so can comfortably sort it here for subsequent binning 

    % define bins
    switch binmethod
        case 'abs'
            [bincounts, binedges2] = histcounts(deps_all,nbins); 
        case 'rank'
            [bincounts, binedges2] = histcounts(1:length(deps_all),nbins, 'BinLimits',[1,length(deps_all)]); 
            binedges2 = deps_all(round(binedges2));
    end
    binedges = nan(nbins,2);
    for k = 1:nbins
        binedges(k,1) = binedges2(k);
        binedges(k,2) = binedges2(k+1);
    end

    % plot bins on histogram
    figure; histogram(deps_all); hold on; 
    for k = 1:nbins+1
        plot( [binedges2(k),binedges2(k)], ylim ,'b--')
    end
    title([ sub ' ' binmetric ' n = ' num2str(length(deps_all)) ' (' sprintf('%d, ',bincounts)  ')' ]);
    
    
    %% loop through file types & load selected data 
    for ft = 1:length(file_types)
        file_type = file_types{ft};
    
        %% get session filenames
        cd(homedir)
        if strcmp(file_types{ft},'EEG') || strcmp(file_types{ft},'GAM')
            
            files = dir; files = {files.name}; 
            files = files( endsWith(files,[ sub '.set']) & ...
                ~contains(files,'all') & ~contains(files,'bin') )';
            files = files(  contains(files,['merged_' file_types{ft}]) );

            % check number of files match the session list
            if length(files) ~= nfiles
                error 'MISMATCH IN FILE-NUMBERS, MISSING SESSIONS OR NEED TO HANDLE FILE-FILTERING?';
            end

            % check that these are the correct sessions
            temp = cell2mat(regexp( files, 'S\d\d' ));
            if length(unique(temp)) ~= 1, error 'files with same format should have session code in the same place...'; end
            temp2 = cellfun( @(x) x(temp(1):temp(1)+2) ,files ,'Uniformoutput',false);
            if ~isequal(temp2,sesh)
                error 'file session codes don''t match session list'
            end

            files2load = files;

        else
            files2load = repmat({['merged_' file_type ' S##_MT ' sub '.set']}, length(sesh), 1);
            files2load = cellfun( @(x,y) strrep(x,'S##',y) ,files2load, sesh, 'UniformOutput', false);
        end
        
        %% load sessions for this file_type
        for f = 1:nfiles
            if strcmp(file_type,'MUA') || strcmp(file_type,'LFP') || strcmp(file_type,'GAM')
                data.(file_type)(f) = pop_loadset('filename',files2load{f},'filepath',s.savePath ...
                    , 'loadmode', 1:5); % only load actual electrodes, not avg of electrodes
            else
                data.(file_type)(f) = pop_loadset('filename',files2load{f},'filepath',s.savePath);
            end
        end
        fprintf('%s ... finished loading\n', file_type)
        
    %% END loop through file-types
    end

    %% get depths of each loaded electrode
    inds = cellfun( @(x) find(strcmp(x, tbl_depths.sesh))  , sesh );
    deps = [];
    for k = 1:5
        temp = tbl_depths.([binmetric '_' sub(end) num2str(k)]);
        deps = [ deps, temp(inds) ];
    end

    %% bin depths
    depbins = arrayfun(@(x) find( x > binedges(:,1) & x <= binedges(:,2)  )  , deps, 'UniformOutput'  , false);
    depbins(cellfun(@isempty,depbins)) = {nan};
    depbins = cellfun(@(x) x , depbins);

    %% loop through file-types
    for ft = 1:length(file_types)
        file_type = file_types{ft};

        %% bin data
        tempout_mean    = cell(nbins,length(s.conds));
        tempout_trials  = cell(nbins,length(s.conds));

        % loop through bins
        for b = 1:nbins
            
            % loop through sessions and assign electrodes to bins
            for f = 1:nfiles
    
                % find & extract any electrodes belonging to this bin in this file
                chans2take = find(depbins(f,:) == b);
                if ~isempty(chans2take) 
                    evalc('temp = pop_select( data.(file_type)(f), ''channel'', chans2take)');  % ? NOTE that this assumes we are dealing with LFP or MUA (will need something else to do EEG, EYE, or MISC)
                    
                    % get each condition
                    for cond = 1:3

                        % extract this condition
                        evalc 'tempcond = pop_select(temp, ''trial'', find(strcmp({temp.event.type},s.conds{cond})) );';
                        tempcond = tempcond.data;
                        tempcond_r = reshape(permute(tempcond,[3 1 2]), [size(tempcond,3)*size(tempcond,1), size(tempcond,2)]);
    
                        tempout_mean{b,cond}    = [tempout_mean{b,cond};   mean(tempcond,3) ];
                        tempout_trials{b,cond}  = [tempout_trials{b,cond}; tempcond_r      ];

                    end % condition loop
                end
    
            end % file loop
            b
        end % bin loop


        %% import binned data into EEGLAB 
        tempd = data.(file_type)(1);
        nsamps = size(tempout_mean{b,cond},2);
        for cond = 1:3
            for b = 1:3
                
                data_binned_means(b,cond).(file_type)   = ...
                    pop_importdata('data', permute(shiftdim(tempout_mean{b,cond},-1),[1 3 2]) ...
                    ,'nbchan', 1, 'pnts',nsamps, 'srate', tempd.srate  ...
                    ,'xmin', tempd.xmin);
                data_binned_trials(b,cond).(file_type) = ...
                    pop_importdata('data',permute(shiftdim(tempout_trials{b,cond},-1),[1 3 2]) ...
                    ,'nbchan', 1, 'pnts',nsamps, 'srate', tempd.srate  ...
                    ,'xmin', tempd.xmin);
                
            end
        end
    
        %% SAVE in EEGLAB
        for cond = 1:3
            for b = 1:3
                
                % save means
                filename_means{cond,b} = [ 'chanavgs' ' ep_' s.conds{cond} ' bin_' num2str(b)  ...
                    ' bm_' binmethod ' '  binmetric ' merged_' file_type ' ' sub  ];
                pop_saveset( data_binned_means(b,cond).(file_type), ...
                    'filepath',tempd.filepath, 'filename', filename_means{cond,b});

                % save trials
                filename_trials{cond,b} = [ 'alltrials' ' ep_' s.conds{cond} ' bin_' num2str(b)  ...
                    ' bm_' binmethod ' '  binmetric ' merged_' file_type ' ' sub  ];
                pop_saveset( data_binned_trials(b,cond).(file_type), ...
                    'filepath',tempd.filepath, 'filename', filename_trials{cond,b});
                
                
            end
        end
    
        %% EXPORT to letswave
        cd([tempd.filepath filesep 'lw'])
        mets = {'means','trials'};
        for cond = 1:3
            for b = 1:3

                for m = 1:2
                    met = mets{m};

                    % get dataset
                    eval(['EEG = data_binned_' met '(b,cond).(file_type);' ])

                    % downsample for lighter lw storage
                    EEG = pop_rs_downsample(EEG,s.dsf);
                
                    % crop for lw
                    EEG = pop_select(EEG, 'time', s.xlims.lw.(file_type) );

                    % rename channel to CZ for superimposing
                    if ~strcmp(file_type,'EEG')
                        EEG.chanlocs(1).labels = 'CZ';
                    end
                
                    % export to lw
                    filename = eval(['filename_' met '{cond,b}']);
                    if strcmp(file_type,'EEG')
                        rs_convert_lab2lw_V1( EEG , filename, chanlocs_lw );
                    else
                        rs_convert_lab2lw_V1( EEG , filename, [] );
                    end 

                end
                
            end
        end


    %% END file type loop
    end 

    
    


%% END loop through subs
end

%%
cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

