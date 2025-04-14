% 
% 
%       plot averages of EEG using LME estimates and confidence intervals 
%           - this version also includes a summary metric of the 
%               LFP pooled electrode average as a predictor to demonstrate relationship between pupil response & VP
%           
%           !!! missing topography plots from LME_EEG
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
split_range  = {[0, 0.2],[0, 0.2],[0, 0.4]}; % range to compute ERP amplitude within

s.conds = {'AUD','SOM','VIS'}; nconds = length(s.conds);

% subs = {'SubM','SubT'};
% subs = {'SubM'};
subs = {'SubT'};

DSF = 2; % downsample factor, 2 = 512, 4 = 256

%% get topoplot stuff

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
    cfg.zscore_win.LFP = s.xlims.plot.LFP;
    cfg.zscore_cond   = true;

    [data,tbl_depths] = VPmonkey_mergeSesh(cfg);

    % downsample for speed
    data.LFP = pop_rs_downsample(data.LFP,DSF);
    data.EEG = pop_rs_downsample(data.EEG,DSF);

    % extract time-window of interest
    data.LFP = pop_select(data.LFP, 'time', s.xlims.plot.LFP);
    data.EEG = pop_select(data.EEG, 'time', s.xlims.plot.EEG);
    nsesh = length(data.LFP);

    %% LOOP THROUGH CONDITIONS
    for cond = 1:length(s.conds)

        %% EXTRACT DATA CORRESPONDING TO CONDITION
       
        %
        clear dcond
        for sh = 1:nsesh
            ents_temp = data.LFP(sh).event;
            inds = find(strcmp({ents_temp.type},s.conds{cond}));
            if ~isempty(inds)
                dcond.LFP(sh) = pop_select(data.LFP(sh), 'trial', inds );
                dcond.EEG(sh) = pop_select(data.EEG(sh), 'trial', inds );
            end
            ntrials_per_sesh(sh) = length(inds);
        end

        % extract all trials
        plotchan = find(strcmpi({data.EEG(1).chanlocs.labels},'Cz'));
        deeg = cat(3,dcond.EEG.data);
        deeg = squeeze(deeg(plotchan,:,:))';
        dlfp = squeeze(cat(3,dcond.LFP.data))';
        ents = cat(2,dcond.LFP.event);
        if length(ents) ~= size(deeg,1)
            error 'mismatch between events and ntrials - need to exclude repeat events'
        end
        [ntrials,nsamps] = size(deeg);
        times = dcond.EEG(1).times/1000;

        %% GET INFO FOR BUILDING PREDICTOR MATRICES

        % get per-trial session list 
        sesh = [];
        for sh = 1:nsesh
            sesh = [sesh; repmat( sh, ntrials_per_sesh(sh), 1) ];
        end 

        % get per trial LFP amplitude
        times_lfp = dcond.LFP(1).times/1000;
        inds = findnearest(times_lfp,split_range{cond}(1)):findnearest(times_lfp,split_range{cond}(2));
        lfp = double(mean(abs(dlfp(:,inds)),2));
        lfpz = (lfp-mean(lfp))/std(lfp,1);

        % binned-LFP amplitude
%         nbins = 3;
%         [~,temp] = sort(lfpz); % applied to rank values but still doesn't ensure equal bin sizes... especially for nbins = 2
%         [lfpz_bin,lfpz_bin_edges] = discretize(temp,nbins,'categorical'); 
        med = median(lfpz); nbins = 2;
        lfpz_bin = dummyvar(categorical(lfpz > med));

        % LFP rank
        [~,lfp_rank] = sort(lfp);
        lfp_rankz = (lfp_rank-mean(lfp_rank))/std(lfp_rank,1);

        % random effects
        G = categorical(sesh);
        Z = ones(ntrials,1);

        %% PREPARE SET OF MODELS

        % build LME predictor matrices
        models = [];
        %
        models(end+1).name    = 'INT'; %#ok<*SAGROW> 
        models(end).x     = ones(ntrials,1);
        models(end).xvars = {'Intercept'};
        models(end).formula = 'EEG ~ 1 + (1 | s)'; %
        % 
        models(end+1).name    = 'LFPz_bin';
        models(end).x     = [ lfpz_bin ]; %#ok<*NBRAK2> 
        models(end).xvars = arrayfun(@(x) ['LFPz_bin_' num2str(x)] , 1:nbins, 'UniformOutput',false);
        models(end).formula = 'EEG ~ lfpz_bin + (1 | s)'; %
        %

        nmodels = length(models);
        for m = 1:nmodels
            models(m).nx = length(models(m).xvars);
        end

        %% LOOP THROUGH MODELS
        for m = 1:nmodels

            X     = models(m).x;
            xvars = models(m).xvars;
            nx    = models(m).nx;

            % PREPARE MODEL OUTPUTS
            est = nan(nx,nsamps);
            upp = nan(nx,nsamps);
            bic = nan(1,nsamps);
    
            %% LOOP THROUGH TIMEPOINTS
            fprintf('%s',repmat('.',1,nsamps))
            fprintf('\n\n')
            tin = tic;
            parfor t = 1:nsamps
                
                %% FIT MODEL
                amp = double(deeg(:,t)); 
   
                % fit linear mixed-effects model
                mdl = fitlmematrix( X, amp, Z, G, ...
                    'FitMethod','REML', ...
                    'FixedEffectPredictors',xvars);
                
                % store outputs
                est(:,t)   = mdl.Coefficients.Estimate;
                upp(:,t)   = mdl.Coefficients.Upper;
                bic(:,t)   = mdl.ModelCriterion.BIC;
            
                fprintf('\b|\n');
%                 fprintf('|'); % ? doesn't work in parfor loop
            % TIME LOOP
            end
            fprintf('\nLME fitting complete (%d/%d) in %.2f mins\n\n',m,nmodels,toc(tin)/60)

            % STORE OUTPUTS
            models(m).est = est;
            models(m).bic = bic;
            models(m).moe = upp - est; % COMPUTE MARGIN OF ERROR (RELATIVE CONFIDENCE INTERVAL)

        %% MODEL LOOP
        end

        %% SAVE MODEL RESULTS
        resfile = [resdir filesep sub '_results_LME_EEG_LFP_' s.conds{cond} '.mat'];
        save(resfile,'models')

        %% (OPTIONAL) LOAD MODEL RESULTS
        if ~exist('models','var')
            load(resfile)
        end

        %% PLOT - MODEL BIC FITS
        figname = [ sub '_LME_EEG_LFP_BIC_' s.conds{cond} ];

        % compute relative BIC values for easier visualisation
        bic_all = cat(1,models.bic);
%         bic_all = bic_all - mean(bic_all); 
        bic_all = bic_all - min(bic_all); % ? this is more interpretable than relative to mean

        % plot
        fig = figure('name',figname); 
        subplot(2,1,1);
        plot(times, mean(deeg),'LineWidth',2);
        subplot(2,1,2);
        plot(times, bic_all, 'LineWidth',2);
        legend(strrep({models.name},'_','-'));
        %
        rs_saveExportFig(fig, figsdir, figname, true);

        %% LOOP THROUGH MODELS AND PLOT ALL COEFFICIENTS
        for m = 1:nmodels
            mdl = models(m);
            %
            figname = [ sub '_LME_EEG_LFP_COEF_' s.conds{cond} '_mdl-' mdl.name  ];

            est2plot = mdl.est;
            moe2plot = mdl.moe;

            % plot coefficient estimate w/ confidence intervals
            fig = figure('name',figname); 
            cols = distinguishable_colors(mdl.nx);
            for k = 1:size(est2plot,1)
                boundedline(times, est2plot(k,:), moe2plot(k,:),...
                    'cmap',cols(k,:),'LineWidth',1); hold on;
            end
            set_boundedline_transparency(0.8);
            xlim(s.xlims.plot.EEG);
            ylim(s.ylims.plot.EEG); 
            set(gca,'YDir','reverse')
            plot(xlim,[0,0],'k-')
            xlabel 'time (s)'
            ylabel 'coefficient estimate (A.U.)'

            % legend
            legs = [repmat({''},1,mdl.nx), strrep(mdl.xvars,'_',' ')];
            [~,inds] = sort([1:mdl.nx 1:mdl.nx]);
            legs = legs(inds);
            legend(legs);
            
            %
            rs_saveExportFig(fig, figsdir, figname, false);
      
        end

        %%
        close all 
        

        %% CONDITION LOOP
    end 

%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
