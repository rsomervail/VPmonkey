% 
% 
%   LME to predict pooled-electrode LFP response with set of EEG IC predictors
% 
%       ** IN THIS VERSION I'm predicting across time, with timepoints as different samples **
% 
%   ! UNFINISHED
%  
%   ! starting to get the feeling that the scores slope coefficients can't be summed the way
%     I want them to... maybe what I actually need is intercepts scaling the mean score?
%           but maybe that would not leverage the actual variability properly?
%       - maybe I should run it across all time-points...? maybe in a restricted window though..
% 
%   ! there are some bad chans on EEG trials to be removed - see  figure; plot(sqz(deeg(1,:,:))) for subM
%       maybe trial-wise interpolation of anything bigger than 200 or 300 uV?
% 
%   ! try decreasing stopping criterion and/or learning rate and seeing if this affects results
%       - check also IC topos, maybe they look better
% 
%   ! right now only running one model with all predictors
%     need a better framework for this analysis, since predictors vary by timepoint now too..
%       - also assuming all predictors when computing the back-projected signal via the scores
% 
%       - Richard Somervail, 2025
%%
clc
clearvars
close all

homedir = '/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_main';
cd(homedir);

figsdir =  [ getRoot '/VPmonkey/paper/figures/raw/' ]; 
resdir =  [ getRoot '/VPmonkey/paper/results/' ]; 

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

%% SETTINGS
s = [];
s.savePath =  [ getRoot '/VPmonkey/paper/results/lw' ];  % mkdir(s.savePath)

% plot limits  
lim = VPmonkey_fetchLimits;

% split_range  = [0, 0.35]; % range to compute ERP amplitude within
split_range  = {[0, 0.3],[0, 0.3],[0, 0.4]}; % range to compute ERP amplitude within

s.conds = {'AUD','SOM','VIS'}; nconds = length(s.conds);

subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};

DSF = 2; % downsample factor, 2 = 512, 4 = 256

%% get topoplot stuff

topo = VPmonkey_fetchTopoPeaks;

s.gridscale = 200; % 70 for quick topoplots, 200 fine for figure-ready plots

addpath([  getRoot filesep 'MATLAB' filesep 'cbrewer' ])

[colormap2]=cbrewer('div', 'RdBu', 1000, 'linear'); % cubic , linear
colormap2 = colormap2(length(colormap2):-1:1, : );

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    savestr = []; % tracking which modules are used to reject trials

    %% LOAD DATA

    % LFP electrode rejection on the basis of unusually high SD value
    ar_lfp_SD.thresh = 1.5;

    % by electrode (LFP, MUA etc)
    cfg = struct; cfg.exportlw = false;
    cfg.sub = sub;
    if exist('ar_lfp_SD','var'), cfg.ar_lfp_SD = ar_lfp_SD; end
    cfg.filetypes = {'LFP','EEG'};
    cfg.LFP_include = {'LP_30'};
%     cfg.LFP_exclude = {'LP_30'};
    cfg.filt_preproc = 'allSesh'; % 'allGood' or 'allSesh'  % ? better with all sessions because more data
    cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
    cfg.EEG_exclude = {'noICA'};
    cfg.byElectrode = false;
    cfg.average   = false;
    cfg.mergesesh = false;

    % z-score
    cfg.zscore = {'LFP'};  % ? here z-scoring seems helpful for improving interpretation of coefficients (also because overall LFP amplitude varies with electrode location/depth)
    cfg.zscore_win.LFP = lim.xlims.plot.LFP;
    cfg.zscore_cond   = true;

    [data,tbl_depths] = VPmonkey_mergeSesh(cfg);

    % downsample for speed
    data.LFP = pop_rs_downsample(data.LFP,DSF);
    data.EEG = pop_rs_downsample(data.EEG,DSF);

    % extract time-window of interest
    data.LFP = pop_select(data.LFP, 'time', lim.xlims.plot.LFP);
    data.EEG = pop_select(data.EEG, 'time', lim.xlims.plot.EEG);
    nsesh = length(data.LFP);

    %% COMPUTE ICA ON EEG DATA ACROSS ALL TRIALS/CONDITIONS OF ALL SESSIONS

    % compute within-condition averages per-session
    data.EEG_avg = pop_rs_average(data.EEG, true);

    % format data
    deeg_avg = double(cat(3,data.EEG_avg.data));
    [nchans,nsamps,nepochs_ica] = size(deeg_avg);
    deeg_avg_r = reshape(deeg_avg,nchans,nsamps*nepochs_ica);

    % run ICA
    [weights,sphere,compvars] = runica(deeg_avg_r);
%     [weights,sphere,compvars,~,~,~,scores_all_r] = runica(deeg_r);
    um = weights*sphere; % um is comps x chans
    mm = inv(um); % mm is chans x comps
    ncomps = length(compvars);

    % compute scores across all trials by applying filter to individual trials
    deeg    = double(cat(3,data.EEG.data));
    ntrials_all  = size(deeg,3);
    deeg_r  = reshape(deeg,nchans,nsamps*ntrials_all);
    scores_all_r = um * deeg_r;
    
    % reshape scores and check they match LFP number of trials
    scores_all = reshape(scores_all_r,ncomps,nsamps,ntrials_all);
    ents_all = cat(2,data.LFP.event); % compute all event codes for later
    if size(scores_all,3) ~= length(ents_all)
        error 'scores & event ntrials mismatch'
    end

    %% PLOT ALL COMPONENT TOPOGRAPHIES
    for c = 1:ncomps
        subplot(6,5,c);
        vals = mm(:,c);
        ctemp = [];
        ctemp.clims     = 'absmax';
        ctemp.colormap  = colormap2;
        ctemp.lay       = data.EEG(1).chanlocs;
        ctemp.gridscale = s.gridscale; % 700
        ctemp.numcontours = 0;
        rickoplot_MK( vals, ctemp );
    end
    drawnow

    %!!! SAVE THIS FIGURE

    %% LOOP THROUGH CONDITIONS
    for cond = 1:length(s.conds)

        %% EXTRACT DATA CORRESPONDING TO CONDITION
       
        % select trials for this condition
        clear dcond
        for sh = 1:nsesh
            ents_temp = data.LFP(sh).event;
            inds = find(strcmp({ents_temp.type},s.conds{cond}));
            if ~isempty(inds)
                dcond.LFP(sh) = pop_select(data.LFP(sh), 'trial', inds );
            end
            ntrials_per_sesh(sh) = length(inds);
        end

        % extract all trials
        dlfp = squeeze(cat(3,dcond.LFP.data))';
        ents = cat(2,dcond.LFP.event);
        if length(ents) ~= size(dlfp,1)
            error 'mismatch between events and ntrials - need to exclude repeat events'
        end
        [ntrials,nsamps] = size(dlfp);
        dlfp_r = double(reshape(dlfp,ntrials*nsamps,1));
        times = dcond.LFP(1).times/1000;

        % extract corresponding trials from scores
        inds = strcmp({ents_all.type},s.conds{cond});
        scores = scores_all(:,:,inds);
        scores = permute(scores,[3,2,1]);
        scores_r = double(reshape(scores,ntrials*nsamps,ncomps));

        %% GET INFO FOR BUILDING PREDICTOR MATRICES

        % get per-trial/sample session list 
        sesh = [];
        for sh = 1:nsesh
            sesh = [sesh; repmat( sh, nsamps*ntrials_per_sesh(sh), 1) ];
        end 

        % random effects
        G = categorical(sesh);
        Z = ones(size(sesh));

        %% PREPARE SET OF MODELS

        % build LME predictor matrices
        models = [];
        %
%         models(end+1).name    = 'INT'; %#ok<*SAGROW> 
%         models(end).x     = ones(ntrials,1);
%         models(end).xvars = {'Intercept'};
%         models(end).formula = 'EEG ~ 1 + (1 | s)'; %
        % 
        models(end+1).name    = 'COMPS_ALL';
        models(end).x     = [ nan(ntrials,ncomps) ]; %#ok<*NBRAK2> 
        models(end).xvars = arrayfun(@(x) ['COMP_' num2str(x)] , 1:ncomps, 'UniformOutput',false);
        models(end).formula = 'EEG ~ comp_1 + comp_2 + ... comp_N + (1 | s)'; %
        %

        nmodels = length(models);
        for m = 1:nmodels
            models(m).nx = length(models(m).xvars);
        end

        %% LOOP THROUGH MODELS
        for m = 1:nmodels

%             X     = models(m).x;
            xvars = models(m).xvars;
            nx    = models(m).nx;

            % PREPARE MODEL OUTPUTS
            est = nan(nx,nsamps);
            upp = nan(nx,nsamps);
            bic = nan(1,nsamps);

            %% FIT MODEL ACROSS ALL TIMEPOINTS AND TRIALS
            tin = tic;
  
            % fit linear mixed-effects model
            mdl = fitlmematrix( scores_r, dlfp_r, Z, G, ...
                'FitMethod','REML', ...
                'FixedEffectPredictors',xvars);
            
            % store outputs
%             est(:,t)   = mdl.Coefficients.Estimate;
%             upp(:,t)   = mdl.Coefficients.Upper;
%             bic(:,t)   = mdl.ModelCriterion.BIC;

            fprintf('\nLME fitting complete (%d/%d) in %.2f mins\n\n',m,nmodels,toc(tin)/60)
    
%             %% LOOP THROUGH TIMEPOINTS
%             fprintf('%s',repmat('.',1,nsamps))
%             fprintf('\n\n')
%             tin = tic;
%             parfor t = 1:nsamps %#ok<PFUIXW> 
%                 
%                 %% FIT MODEL
%                 amp = double(dlfp(:,t)); 
% 
%                 % for all models with comps, replace nans with scores
%                 X = scores(:,:,t);
%    
%                 % fit linear mixed-effects model
%                 mdl = fitlmematrix( X, amp, Z, G, ...
%                     'FitMethod','REML', ...
%                     'FixedEffectPredictors',xvars);
%                 
%                 % store outputs
%                 est(:,t)   = mdl.Coefficients.Estimate;
%                 upp(:,t)   = mdl.Coefficients.Upper;
%                 bic(:,t)   = mdl.ModelCriterion.BIC;
%             
%                 fprintf('\b|\n');
% %                 fprintf('|'); % ? doesn't work in parfor loop
%             % TIME LOOP
%             end
%             fprintf('\nLME fitting complete (%d/%d) in %.2f mins\n\n',m,nmodels,toc(tin)/60)
% 
%             % STORE OUTPUTS
            models(m).est = mdl.Coefficients.Estimate;
            models(m).bic = mdl.ModelCriterion.BIC;
            models(m).moe = mdl.Coefficients.Upper - mdl.Coefficients.Estimate; % COMPUTE MARGIN OF ERROR (RELATIVE CONFIDENCE INTERVAL)

        %% MODEL LOOP
        end

        %% SAVE MODEL RESULTS
        resfile = [resdir filesep sub '_results_LME_LFP_EEG_PNT' s.conds{cond} '.mat'];
        save(resfile,'models')

        %% (OPTIONAL) LOAD MODEL RESULTS
        if ~exist('models','var')
            load(resfile)
        end

        %% WEIGHT SCORES BY ESTIMATED MODEL COEFFICIENTS AND PLOT
        

        
        %%
        close all 
        

        %% CONDITION LOOP
    end 

%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
