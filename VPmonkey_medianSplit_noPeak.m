% 
% 
%       median split of data of N types according to values of another type
%           using summed absolute amplitude within time-window
%               08/11/24 - I think this version was the one I ended up settling on because it was easier to interpret
% 
%
%%
clc
clearvars
close all

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

homedir = ['/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_' subfold];
cd(homedir);

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

pathlw  = [getRoot 'toolbox' filesep 'letswave7-master'];
pathlab = [getRoot 'toolbox' filesep 'eeglab2022.1'];

%% SWITCH TO EEGLAB
evalc 'addpath(genpath(pathlab));';
evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''fieldtrip-20210311'' ]))'; % remove fieldtrip to tidy things up a bit
evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''Biosig3.7.9'' ]))'; % remove biosig to tidy things up a bit
evalc 'rmpath(genpath(pathlw));';

%% SETTINGS
s = [];
s.savePath =  [homedir filesep 'lw']; % mkdir(s.savePath)
s.savePath_figs =  [ getRoot '/VPmonkey/paper/figures/raw' ]; % mkdir(s.savePath_figs)

% epoching 
s.xlims.lw.EEG  = [-0.5 1];  
s.xlims.lw.LFP  = [-0.5 1];  
s.xlims.lw.MUA  = [-0.5 1];  
s.xlims.lw.MISC = [-0.2 0.5];
s.xlims.lw.EYE  = [-0.5 6]; 

% plot limits 
% xlims
s.xlims.plot.EEG  = [-0.2 0.6];  
s.xlims.plot.LFP  = [-0.2 0.6];  
s.xlims.plot.MUA  = [-0.2 0.6];  
s.xlims.plot.EYE  = [-1   3]; 
s.xlims.plot.MISC = [-0.2 0.5];
% ylims
% ylims - amplitudes
s.ylims.plot.EEG  = [-20 20];  
s.ylims.plot.LFP  = [-4 4];   
s.ylims.plot.MUA  = [-0.14 0.14];  
s.ylims.plot.EYE  = [-1 1]; 
% ylims - tvals
s.ylims.tvals.EEG  = [-6 6];  
s.ylims.tvals.LFP  = [-40 40];  
s.ylims.tvals.MUA  = [-6 6]; 
s.ylims.tvals.EYE  = [-6 6];
% topoplot clims - amplitudes
s.clims.SubM.AUD = [-5 5];
s.clims.SubM.SOM = [-5 5];
s.clims.SubM.VIS = [-5 5];
s.clims.SubT.AUD = [-5 5];
s.clims.SubT.SOM = [-5 5];
s.clims.SubT.VIS = [-5 5];

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

file_types   = {'EYE','LFP'}; % never use more than 2 filetypes or there will be clashes with sesh selection
% file_types   = {'MUA','LFP'};
% file_types   = {'EEG','LFP'}; % ? this analysis is really intended for MUA and EYE, not EEG (use CCA for that)

% split_type   = 'EEG'; split_chan   = 'Cz'; 
split_type   = 'LFP'; split_chan = 'LFP'; % ? makes the most sense analytically but haven't tested script out like this yet
% split_range  = [0, 0.3]; % range to compute ERP amplitude within
split_range  = [0, 0.4]; % range to compute ERP amplitude within

s.filt_elec = '';
% s.filt_elec = 'Depth > MED'; % excluding median electrode here so that num electrodes is equal
% s.filt_elec = 'Depth < MED';
% s.filt_elec = 'DepthRelTOA > MED';
% s.filt_elec = 'DepthRelTOA < MED';

%% import intracranial information table
tbl_depths = VPmonkey_import_electrodeDepths;

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};

    %% LOAD SELECTED DATA
    % by session (EEG)
    cfg = struct;
    cfg.sub = sub;
        % trial-rejection criteria
        ar = struct;
        ar.method = 'time';
        ar.metric = 'median';
        ar.thresh = 3;
        ar.timeprop = 0.15; % ?
        ar.chanprop = 0.15; % ?
        cfg.autoar = ar;
    cfg.filetypes = file_types;
    cfg.zscore    = file_types;
    cfg.zscore_cond = true;
    cfg.byElectrode = true;
    cfg.average   = false;
    cfg.mergesesh = true;
    cfg.filt_elec = s.filt_elec;
    if any(contains(file_types,'EEG'))
        cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
        cfg.EEG_exclude = {'noICA'};
        cfg.filt_preproc = 'allGood';  % these refer to EEG preproc
    else
        cfg.filt_preproc = 'allSesh';
    end
    if any(contains(file_types,'EYE'))
        cfg.filt_eye = true; % only get good EYE sessions
        cfg.EYE_include = {'LP_5'}; % only include low-pass filtered pupil data
    end
    [data, tbl_depths] = VPmonkey_mergeSesh(cfg);
    
    fprintf('... finished loading data\n')%

    %% extract data channel to be used for splitting
    split_chan_ind = find(strcmpi(  {data.(split_type).chanlocs.labels}, split_chan   ));
    data.split = pop_select( data.(split_type), 'channel', split_chan_ind);

    %% low-pass filter data used to compute median splits
    data.split = pop_eegfiltnew( data.split, 'hicutoff', 30 );

    %% high-pass filter data used to compute median splits
    data.split = pop_eegfiltnew( data.split, 'locutoff', 1 );

    %% copy filtered version to it's main field for later export of appropriately filtered signal
    data.(split_type) = data.split;

    %% get only certain time-range of data
    data.split = pop_select( data.split, 'time',  split_range);

    %% LOOP THROUGH CONDITIONS
    data.split = pop_rs_removeNonZeroEvents(data.split);
    for cond = 1:length(s.conds)

        %% concatenate data to be used for splitting
        d = data.split.data;
        ents = {data.split.event.type};
        d(:,:,~strcmp(ents,s.conds{cond})) = []; % remove epochs which correspond to different modality
        d = permute(d,[3 2 1]); % trials, samples, channels (? nchans=1)

        %% compute mean absolute signal within window, per trial
        
        % compute average
        dm = mean(d);

        % compute per-trial mean absolute ERP amplitude across window
        pks = mean(abs(d),2); % simple version of mean absolute amplitude
        pkstr = 'meanAbs';

%         % compute per-trial mean absolute ERP amplitude across window
%         % (flips signs of timepoints which are negative, allowing reverse polarity peaks in trials to be considered as "low amplitude")
%         pksigns = ones(size(dm));
%         pksigns(dm<0) = -1;
%         pks = mean( d.*pksigns, 2); 
%         pkstr = 'meanSignFlipped';
            
        %% median split
        [~, inds] = sort( pks );
        mid = floor(length(inds)/2);    

        %% loop through file types and median split
        for ft = 1:length(file_types)
            file_type = file_types{ft};

            % get this condition
            dtemp = data.(file_types{ft});
            dtemp = pop_rs_removeNonZeroEvents(dtemp);
            dtemp = pop_select(dtemp,...
                'trial', find(strcmp({dtemp.event.type},s.conds{cond})));

            % split 
            for m = 1:2 % low to high
                if m == 1
                    dtemp_bin  = pop_select(dtemp,'trial', inds(1:mid) ); 
                    msplitcode = 'M_LOW';
                else
                    dtemp_bin  = pop_select(dtemp,'trial', inds(mid+1:end) ); 
                    msplitcode = 'M_HIGH';
                end

                % EXPORT TO LW
                % downsample for efficiencyfile_type
                if ~strcmp( file_type, 'MISC' )
                    dtemp_bin  = pop_rs_downsample(dtemp_bin,s.dsf);
                end
    
                % crop for lw
                dtemp_bin = pop_select(dtemp_bin, 'time', s.xlims.lw.(file_type) );
    
                % replace all NaNs with zeros for easier lw plotting
                dtemp_bin.data(isnan(dtemp_bin.data)) = 0;
        
                % export to lw
                cd(s.savePath)

                % define electrode filter string
                if ~exist('filt_elec_str','var') && ~isempty(s.filt_elec)
                    filt_elec_str = ['filtElec_YES ' strrep(s.filt_elec,' ','_')];
                    filt_elec_str = strrep(filt_elec_str,'>=','MTEq');
                    filt_elec_str = strrep(filt_elec_str,'<=','LTEq');
                    filt_elec_str = strrep(filt_elec_str,'>','MT');
                    filt_elec_str = strrep(filt_elec_str,'<','LT');
                elseif ~exist('filt_elec_str','var') && isempty(s.filt_elec)
                    filt_elec_str = 'filtElec_NO';
                end
                % build saveName
                saveName = [ 'ep_' s.conds{cond} ...
                    ' medsplit_noPeak_s' split_type '_' strjoin(file_types(~strcmp(file_types,split_type)),'_')  ...
                    ' ' filt_elec_str ...
                    ' ' pkstr  ...
                    ' merged_' file_type ...
                    ' ' msplitcode ...
                     ' ' sub ];
                if strcmp(file_type,'EEG')
                    rs_convert_lab2lw_V1( dtemp_bin ...
                        , saveName, chanlocs_lw);
                else
                    rs_convert_lab2lw_V1( dtemp_bin ...
                        , saveName, []);
                end

            end % median split loop

            %% run cluster-permutation test between low and high splits
            ctemp = struct;
            ctemp.nperms = 100;
            warning 'CURRENTLY ONLY RUNNING 100 PERMS FOR SPEED'
    
            ctemp.alpha1 = 0.05; % could lower this threshold if no clusters are formed, without affecting false positive rate
            ctemp.alpha2 = 0.05; % alpha-level for cluster significance - do not change
    
%                 clocal.cluster_metric = 'maxT';  % fewer false positives
            ctemp.cluster_metric = 'meanT'; % fewer false negatives 
    
            ctemp.minnumchan = 0; % ? for now allowing single-channel clusters and hoping their cluster-T will be low
    %                                  if there are many spurious single-channel clusters, then increase this
    
            if dtemp.nbchan > 1
                ctemp.neighbours = 0.8;  % produces an average of ~9 channels
            end

            % get data
            dtemp_low   = pop_select(dtemp,'trial', inds(1:mid) ); 
            dtemp_high  = pop_select(dtemp,'trial', inds(mid+1:end) ); 

            % RUN CLUSTER PERMUTATION TWO-SAMPLE T-TEST
            [cluster_tvals, cluster_pvals, cluster_clusterids] = pop_rs_clusterperm_ttest2(dtemp_high,dtemp_low,ctemp);

            %% PLOT MEDIAN SPLIT RESULTS

            % get plot channel c
            if strcmp(file_type,'EEG')
                c = find(strcmpi({dtemp.chanlocs.labels},'cz'));
            else
                c = 1; % for all file-types can assume only one channel by this point
            end
            plotchan = dtemp.chanlocs(c).labels;

            % extract data for this channel
            dlow  = dtemp_low.data(c,:,:);
            dhigh = dtemp_high.data(c,:,:);

            % compute means
            dlowm  = mean(dlow,3);
            dhighm = mean(dhigh,3);

            % extract t-value and p-value timecourses for cluster permutation t-test for this channel
            tvals           = squeeze(cluster_tvals.data(c,:,1)); 
            tvals_thresh    = squeeze(cluster_tvals.data(c,:,2)); 
            pvals           = squeeze(cluster_pvals.data(c,:,1));
            pvals_thresh    = squeeze(cluster_pvals.data(c,:,2)); 
            % IF I CHANGE THIS TO PLOT SIG CLUSTERS, NEED TO RUN MORE PERMUTATIONS

            % PLOT - averages
            figname = [ 'plot-medsplit_noPeak-AVG ' sub ' ' s.conds{cond} ' s' split_type '_' file_type ' ' pkstr ]
            figs.AVG = figure('name',figname, 'NumberTitle','off');
            %
            times = dtemp_low.times/1000;
            plot(times, dlowm,  'r', 'LineWidth',2.5 ); hold on;
            plot(times, dhighm, 'b', 'LineWidth',2.5 ); hold on;
            xlim(s.xlims.plot.(file_types{ft}));                    
            ylim(s.ylims.plot.(file_types{ft}));
            plot(xlim,[0,0],'k')

            % PLOT - cluster t plots - thresholded and unthresholded
            figname = [ 'plot-medsplit_noPeak-TVALS ' sub ' ' s.conds{cond} ' s' split_type '_' file_type ' ' pkstr ]
            figs.TVALS = figure('name',figname, 'NumberTitle','off');
            %
            plot(times, tvals, 'k', 'LineWidth', 1); hold on;
            plot(times, tvals_thresh, 'k', 'LineWidth', 2.5); hold on;
            xlim(s.xlims.plot.(file_types{ft}));                    
            ylim(s.ylims.tvals.(file_types{ft}));
            plot(xlim,[0,0],'k')

            % PLOT - cluster p plots - thresholded and unthresholded
            figname = [ 'plot-medsplit_noPeak-PVALS ' sub ' ' s.conds{cond} ' s' split_type '_' file_type ' ' pkstr ]
            figs.PVALS = figure('name',figname, 'NumberTitle','off');
            %
            plot(times, pvals, 'k', 'LineWidth', 1); hold on;
            plot(times, pvals_thresh, 'k', 'LineWidth', 2.5); hold on;
            xlim(s.xlims.plot.(file_types{ft}));                    
            ylim([0,1]);
            plot(xlim,[0.05,0.05],'k')

            % SAVE & EXPORT ALL CURRENT FIGURES
            cd(s.savePath_figs)
            fignames  = fields(figs);
            nfigs = length(fignames);
            for f = 1:nfigs
                rs_saveExportFig(figs.(fignames{f}), s.savePath_figs, figs.(fignames{f}).Name);
            end
            clear figs
            close all

        %% END LOOP THROUGH FILE TYPES
        end 

    %% END LOOP THROUGH CONDITIONS
    end

%% END loop through subs
end

%%
cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

