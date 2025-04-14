% 
% 
%       plot averages of all data types across all sessions
%       using LME estimates and confidence intervals for EEG ONLY
%          ? for LFP I need a separate script which compares to site and depth models using BIC (or WAIC if possible!!)
% 
%  !!! UNFINISHED 
% 
%   - TRY WITH RANDOM EFFECTS it might be better to use lme anyway, with random effects for different sessions
%           - in timepoints where they are negligible, what to do? use LM instead?
% 
%       - it seems that overfitted near-zero random effects might artificially reduce
%           the estimate of the model uncertainty, so I need to iteratively
%           design models using BIC to see whether the random effects are needed
% 
%       - Richard Somervail, 2023
%%
clc
clearvars
close all

homedir = '/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_main';
cd(homedir);

addpath([getRoot '/VPmonkey/Giaplotc_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])


%% SETTINGS
s = [];
s.savePath_figs =  [ getRoot '/VPmonkey/paper/figures/raw' ]; % mkdir(s.savePath_figs)
s.savePath =  [ getRoot '/VPmonkey/paper/results/lw' ];  % mkdir(s.savePath)

% plot limits 
% xlims
s.xlims.plot.EEG  = [-0.2 0.6];  
s.xlims.plot.LFP  = [-0.2 0.6];  
s.xlims.plot.MUA  = [-0.2 0.6];  
s.xlims.plot.EYE  = [-1   3]; 
s.xlims.plot.MISC = [-0.2 0.5];
% ylims - amplitudes
s.ylims.plot.EEG  = [-20 20];  
s.ylims.plot.LFP  = [-100 100];   
% s.ylims.plot.LFP  = [-4 4];   
s.ylims.plot.MUA  = [-0.14 0.14];  
s.ylims.plot.EYE  = [-2e-2 2e-2]; 
% ylims - tvals
s.ylims.tvals.EEG  = [-40 40];  
s.ylims.tvals.LFP  = [-80 80];  
s.ylims.tvals.MUA  = [-40 40]; 
s.ylims.tvals.EYE  = [-40 40];
% topoplot clims - amplitudes
s.clims.SubM.AUD = [-5 5];
s.clims.SubM.SOM = [-5 5];
s.clims.SubM.VIS = [-5 5];
s.clims.SubT.AUD = [-5 5];
s.clims.SubT.SOM = [-5 5];
s.clims.SubT.VIS = [-5 5];

% channels to plot
s.chans2plot.EEG = {'CZ'};
s.chans2plot.MUA = {'MUA'};
s.chans2plot.LFP = {'LFP'};
s.chans2plot.EYE = {'eye_PupilDiameter'};

% EEG topo peaks  ? this will find the closest peak on the average to the specified points
s.topopeaks.SubM.AUD = [30 95 160]/1000; 
s.topopeaks.SubM.SOM = [36 97 177] / 1000
s.topopeaks.SubM.VIS = [22 43  67 98 120 210 290 ]/1000;
s.topopeaks.SubT.AUD = [33 79 150]/1000;
s.topopeaks.SubT.SOM = [36 120 170] / 1000
s.topopeaks.SubT.VIS = [19 46  110  190 290 ]/1000;

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

s.conds = {'AUD','SOM','VIS'};

subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};

file_types = {'EEG'};
% file_types = {'EYE','EEG'};


%% get topoplot stuff
addpath([  getRoot filesep 'MATLAB' filesep 'cbrewer' ])

[colormap2]=cbrewer('div', 'RdBu', 1000, 'linear'); % cubic , linear
colormap2 = colormap2(length(colormap2):-1:1, : );

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    savestr = []; % tracking which modules are used to reject trials

    %% LOAD DATA

%     % trial-rejection criteria
%     ar = struct;
%     ar.method = 'time';
%     ar.metric = 'median';
%     ar.thresh = 3;
%     ar.timeprop = 0.1; % ? too strict?
%     ar.chanprop = 0.1; % ? too strict?

    % by session (EEG, EYE)
    cfg = struct;
    cfg.sub = sub;
%     cfg.autoar = ar;
    cfg.byElectrode = false;
    cfg.average   = false;
    cfg.mergesesh = false;
    % EEG
    if ismember(file_types,{'EEG'})
        cfg.filetypes = {'EEG'};
        cfg.filt_preproc = 'allGoodEEG'; % 'allGoodEEG', 'allPerfEEG', 'allSesh'
        cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
        cfg.EEG_exclude = {'noICA'};
        %
        data_EEG = VPmonkey_mergeSesh(cfg);
        data.EEG = data_EEG.EEG
    end
    % EYE
    if ismember(file_types,{'EYE'})
        cfg.filetypes = {'EYE'};
        cfg.filt_eye = true;
        cfg.EYE_include = {'LP_5'};
        %
        data_EYE = VPmonkey_mergeSesh(cfg);
        data.EYE = data_EYE.EYE;
    end

   

    %% loop through file_types 
    for ft = 1:length(file_types)
        ftype = file_types{ft};

        % extract time-window of interest
        data.(ftype) = pop_select(data.(ftype), 'time', s.xlims.plot.(ftype));

        % extract all trials
        d = cat(3,data.(ftype).data);
        ents = cat(2,data.(ftype).event);
        if length(ents) ~= size(d,3)
            error 'mismatch between events and ntrials - need to exclude repeat events'
        end
        [nchans,nsamps,ntrials] = size(d);
        times = data.(ftype)(1).times/1000;
        chans = {data.(ftype)(1).chanlocs.labels};
        seshlist = [];
        nsesh = length(data.(ftype));
        for k = 1:nsesh
            seshlist = [seshlist; repmat(k,data.(ftype)(k).trials,1)];
        end

        %% LOOP THROUGH CONDITIONS
        conds = unique({ents.type});
        nconds = length(conds);
        for cond = 1:length(s.conds)

            %% EXTRACT DATA CORRESPONDING TO CONDITION
            indscond = strcmp({ents.type}, conds{cond});
            seshlist_cond = seshlist(indscond);
            seshlist_cond_cat = categorical(seshlist_cond);
            dcond    =  d(:,:,indscond);
            entscond = ents(indscond);
            ntrials_cond = sum(indscond);
            ntrials_cond_persesh = arrayfun(@(x) sum(seshlist_cond==x),unique(seshlist_cond));

            %% PREPARE MODEL OUTPUTS
            est = nan(nchans,nsamps,nsesh);
            se  = nan(nchans,nsamps,nsesh);

            % compute critical t-statistic (for confidence intervals later)
            tcrit = tinv(1 - 0.05/2, ntrials_cond - 1);  % t-distribution critical value
      
            %% LOOP THROUGH CHANNELS
            for c = 1:nchans
    %             c = find(strcmpi(chans,s.chans2plot.(ftype){1}));
                dchan = squeeze(dcond(c,:,:));
    
                %% LOOP THROUGH TIMEPOINTS
                for t = 1:nsamps
    %             t = 304; 
                    amp = double(dchan(t,:))';
           
                    %% FIT MODEL

%                     % fit linear model
%                     m = fitlm(seshlist_cond_cat, amp, ...
%                         'intercept',false,'DummyVarCoding','full','RobustOpts',false); % ? iteration limit reached with robustopts, how to increase this? maybe general optimisation settings 
%                     %
%                     est(c,t,:) = m.Coefficients.Estimate;
%                     se(c,t,:)  = m.Coefficients.SE;

                    % fit linear mixed-effects model
                    m = fitlmematrix(seshlist_cond_cat, amp, ...
                        'intercept',false,'DummyVarCoding','full','RobustOpts',false); % ? iteration limit reached with robustopts, how to increase this? maybe general optimisation settings 
                    %
                    est(c,t,:) = m.Coefficients.Estimate;
                    se(c,t,:)  = m.Coefficients.SE;
                
                %% TIME LOOP
                end

                c
            %% CHAN LOOP
            end

            %% MERGE ESTIMATES ACROSS SESSIONS

%             % unweighted
%             est_mean = mean(est,3);
%              ci_mean = mean(ci,3); 
%              % ? should I be weighting these according to number of trials per sesh or
%              %    instead adjust for that weighting so they are balanced...?? probs the latter
% 
%             % weighted by N trials per sample
%             est_mean = nan(nchans,nsamps);
%              ci_mean = nan(nchans,nsamps);
%             for c = 1:nchans
%                 for t = 1:nsamps
%                     est_mean(c,t) = sum(squeeze(est(c,t,:)) .* ntrials_cond_persesh) / sum(ntrials_cond_persesh);
%                     ci_mean(c,t)  = sum(squeeze( ci(c,t,:)) .* ntrials_cond_persesh) / sum(ntrials_cond_persesh);
%                 end
%             end

            % weighted average coefficient & overall margin of error (confidence intervals)
            est_mean = nan(nchans,nsamps);
            moa_mean = nan(nchans,nsamps);
            moa      = nan(size(est)); % separate margin of error for each estimate
            for c = 1:nchans
                for t = 1:nsamps

                    % per estimator
                    moa(c,t,:) = squeeze(se(c,t,:)) .* tcrit;

                    % across estimators
                    w = arrayfun(@(x) 1/(x^2), squeeze(se(c,t,:)));
                    est_mean(c,t) = sum(w .* squeeze(est(c,t,:))) / sum(w);
                    se_comb = sqrt(1/sum(w));
                    moa_mean(c,t) = se_comb * tcrit; % margin of error (combined SE x critical t-value)
                end
            end


            %% PLOT 
            % ? for plot channel only, + topographies of estimate for EEG
            plotchan = strcmpi(chans,s.chans2plot.(ftype){1});
            
            % plot mean coefficient estimate w/ confidence intervals
            figure; 
            boundedline(times, est_mean(plotchan,:), moa_mean(plotchan,:) );
            hold on;
            xlim(s.xlims.plot.(ftype));
            ylim(s.ylims.plot.(ftype));
            plot(xlim,[0,0],'k-')
            xlabel 'time (s)'
            ylabel 'coefficient estimate (uV)'
            %
            figname = [ sub '_LME_avgCI_'  ftype '_' s.chans2plot.(ftype){1} '_' s.conds{cond} ];
            rs_saveExportFig(figname, figsdir, figname);


            %% CONDITION LOOP
        end 

    end
    %% END loop through file_types

%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
