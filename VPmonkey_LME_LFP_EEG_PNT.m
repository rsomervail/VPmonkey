% 
% 
%   LME to predict pooled-electrode LFP response with set of EEG IC predictors
% 
%       ** IN THIS VERSION I'm predicting across time, with timepoints as different samples **
% 
%   ! UNFINISHED - see todo list
% 
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

twin  = [0, 0.4]; % range to compute ERP amplitude within
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
    %
    % lfp only
%     cfg.zscore = {'LFP'}; 
%     cfg.zscore_win.LFP = lim.xlims.plot.LFP;
    % lfp and eeg
    cfg.zscore = {'LFP','EEG'};  
    cfg.zscore_win.EEG = lim.xlims.plot.EEG;
    cfg.zscore_win.LFP = lim.xlims.plot.LFP;
    %
    cfg.zscore_cond    = true;

    [data,tbl_depths] = VPmonkey_mergeSesh(cfg);

    % downsample for speed
    data.LFP = pop_rs_downsample(data.LFP,DSF);
    data.EEG = pop_rs_downsample(data.EEG,DSF);

    % extract time-window of interest
    data.LFP = pop_select(data.LFP, 'time', twin);
    data.EEG = pop_select(data.EEG, 'time', twin);
    nsesh = length(data.LFP);
    times = data.EEG(1).times/1000;

    %% COMPUTE ICA ON EEG DATA ACROSS ALL TRIALS/CONDITIONS OF ALL SESSIONS

    % compute within-condition averages per-session
    data.EEG_avg = pop_rs_average(data.EEG, true);

    % format data
    deeg_avg = double(cat(3,data.EEG_avg.data));
    [nchans,nsamps,nepochs_ica] = size(deeg_avg);
    deeg_avg_r = reshape(deeg_avg,nchans,nsamps*nepochs_ica);

    % run ICA
%     [weights,sphere,compvars] = runica(deeg_avg_r);
    [weights,sphere,compvars] = runica(deeg_avg_r,'extended',1, 'stop', 1e-7); % default for <33 chans is 1e-6
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

    %% COMPUTE COMPONENT METRICS
    % ? not immediately clear whether to use the scores or back-projected EEG GFP
    % ? also not clear whether to use the average response or trial-by-trial (might 
    %    clarify the more selective components, which are likely to be non-phaselocked or artifactual)
    %       - trying average first

    % compute mean absolute score in each condition
    scores_cond = nan(ncomps,nconds);
    for cond = 1:nconds
        trls = find(strcmp({ents_all.type},s.conds{cond}));
        inds = findnearest(times,split_range{cond}(1)):findnearest(times,split_range{cond}(2));
        scores_cond(:,cond) = mean(abs(scores_all(:,inds,trls)),2:3);
    end

    % compute supramodality score
    comp_sup = nan(ncomps,1);
    for c = 1:ncomps
        comp_sup(c) = supramodality_score(scores_cond(c,:));
    end
    [~, ~, comp_sup_rank] = unique(comp_sup); % rank version

    %% PLOT ALL COMPONENT TOPOGRAPHIES
    figname = [ sub '_LME_LFP_EEG_PNT_TOPO-ALLCOMPS' ]; 
    fig = figure('name',figname, 'NumberTitle','off');
    for c = 1:ncomps
        subplot(6,5,c);
        vals = mm(:,c);
        ctemp = [];
        ctemp.clims     = 'absmax';
        ctemp.colormap  = colormap2;
        ctemp.lay       = data.EEG(1).chanlocs;
        ctemp.gridscale = 100; 
        ctemp.numcontours = 0;
        rickoplot_MK( vals, ctemp );
    end

    rs_saveExportFig(fig, figsdir, figname);
    close(fig)

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



        %% FIT MODEL ACROSS ALL TIMEPOINTS AND TRIALS
        
        %!! make sure to store the xvar names

        tin = tic;
        models = [];

        % fit linear mixed-effects model with all predictors
        mdl = fitlmematrix( scores_r, dlfp_r, Z, G, ...
            'FitMethod','REML', ...
            'FixedEffectPredictors',xvars);
        models(end+1).name = 'all_comps';
        models(end).est = mdl.Coefficients.Estimate;
        models(end).bic = mdl.ModelCriterion.BIC;
        models(end).moe = mdl.Coefficients.Upper - mdl.Coefficients.Estimate; % COMPUTE MARGIN OF ERROR (RELATIVE CONFIDENCE INTERVAL)
        fprintf('\nLME fitting complete - %s (%d/%d) in %.2f mins\n\n',models(end).name,m,nmodels,toc(tin)/60)

        % fit a range of models with single component predictors
        % !!! FIND BEST OF ALL OF THESE
        tin = tic;
%         [comp_sup_sorted,sortinds] = sort(comp_sup);
        [comp_sup_sorted,sortinds] = sort(comp_sup,'descend');
        %
        for k = 1:ncomps
            mdl = fitlmematrix( scores_r(:,k), dlfp_r, Z, G, ...
                'FitMethod','REML');
            %       
            models(end+1).name = sprintf('comp_%d',k);
            models(end).est = mdl.Coefficients.Estimate;
            models(end).bic = mdl.ModelCriterion.BIC;
            models(end).moe = mdl.Coefficients.Upper - mdl.Coefficients.Estimate; % COMPUTE MARGIN OF ERROR (RELATIVE CONFIDENCE INTERVAL)
            fprintf('\nLME fitting complete - %s (%d/%d) in %.2f mins\n\n',models(end).name,m,nmodels,toc(tin)/60)
        end   

        % fit a range of models with different supramodality thresholds
        % !!! FIND BEST OF ALL OF THESE
        % ? these models decrease in BIC in either direction, but the rate is different...
%         tin = tic;
% %         [comp_sup_sorted,sortinds] = sort(comp_sup);
%         [comp_sup_sorted,sortinds] = sort(comp_sup,'descend');
%         %
%         for k = 1:ncomps
%             mdl = fitlmematrix( scores_r(:,sortinds(1:k)), dlfp_r, Z, G, ...
%                 'FitMethod','REML');
%             %       
%             models(end+1).name = sprintf('comp_sup_%d',k);
%             models(end).est = mdl.Coefficients.Estimate;
%             models(end).bic = mdl.ModelCriterion.BIC;
%             models(end).moe = mdl.Coefficients.Upper - mdl.Coefficients.Estimate; % COMPUTE MARGIN OF ERROR (RELATIVE CONFIDENCE INTERVAL)
%             fprintf('\nLME fitting complete - %s (%d/%d) in %.2f mins\n\n',models(end).name,m,nmodels,toc(tin)/60)
%         end
        
%         % fit linear mixed-effects models w/ greedy forward model selection
%         % ? make sure to read up on the theory behind this algorithm
%         % ? maybe try an algorithm which prefers fewer predictors?
%         tin = tic;
%         [mdl, ~, selpreds] ...
%             = greedyForwardModelSelection(scores_r, dlfp_r, Z, G);
%         %
%         models(end+1).name = 'greedy_fwd';
%         models(end).est = mdl.Coefficients.Estimate;
%         models(end).bic = mdl.ModelCriterion.BIC;
%         models(end).moe = mdl.Coefficients.Upper - mdl.Coefficients.Estimate; % COMPUTE MARGIN OF ERROR (RELATIVE CONFIDENCE INTERVAL)
%         fprintf('\nLME fitting complete - %s (%d/%d) in %.2f mins\n\n',models(end).name,m,nmodels,toc(tin)/60)

        % count models and compute number of predictors
        nmodels = length(models);
        for m = 1:nmodels
            models(m).nx = length(models(m).xvars);
        end

        %% SAVE MODEL RESULTS
        resfile = [resdir filesep sub '_results_LME_LFP_EEG_PNT' s.conds{cond} '.mat'];
        save(resfile,'models')

        %% (OPTIONAL) LOAD MODEL RESULTS
        if ~exist('models','var')
            load(resfile)
        end

        %% PLOT MODEL TOPOGRAPHIES
        mdl = models(6);

        figname = [ sub '_LME_LFP_EEG_PNT_TOPO-' mdl.name ]; 
        fig = figure('name',figname, 'NumberTitle','off');
%         for c = 1:ncomps  % !!! NEEDS TO LOOP THROUGH ALL COMPS FOR THE mth model
%             subplot(6,5,c);
            vals = mm(:,str2double(mdl.name(end))); %!! do this properly by storing comp indices in models structure
            ctemp = [];
            ctemp.clims     = 'absmax';
            ctemp.colormap  = colormap2;
            ctemp.lay       = data.EEG(1).chanlocs;
            ctemp.gridscale = 100; 
            ctemp.numcontours = 0;
            rickoplot_MK( vals, ctemp );
%         end
    
%         rs_saveExportFig(fig, figsdir, figname);
%         close(fig)

        %% PLOT - BAR CHART OF ALL MODEL BIC VALUES
        %!!!

        %% !!!!!!!!!!! ADAPT BELOW FIGURES FOR NEW MODEL SET

        %% WEIGHT SCORES BY ESTIMATED MODEL COEFFICIENTS AND PLOT
%         for m = 1:nmodels
            m = 1;
            mdl = models(m);
            %
            figname = [ sub '_LME_LFP_EEG_PNT_COEF_' s.conds{cond} '_mdl-' mdl.name  ];

            % weight scores by estimated model coefficients
            scores_w_r = scores_r .* mdl.est';
            
            % reformat to trials and plot means
            scores_w = reshape(scores_w_r,ntrials,nsamps,mdl.nx);
            scores_w_m = squeeze(mean(scores_w,1));

            % plot coefficient estimate scaled mean scores
            fig = figure('name',figname); 
            cols = distinguishable_colors(mdl.nx);
            subplot(2,1,1);
            for c = 1:mdl.nx
                plot(times, scores_w_m(:,c), 'Color', cols(c,:), 'linewidth',2); hold on;
            end
            xlim(lim.xlims.plot.LFP);
%             ylim([-0.2 0.2]); 
            set(gca,'YDir','reverse')
            plot(xlim,[0,0],'k-')
            xlabel 'time (s)'
            ylabel 'amplitude (uV)'
            
            % legend
            % if not plotting CI
            legend(strrep(mdl.xvars,'_',' '));
            % if plotting CI
%             legs = [repmat({''},1,mdl.nx), strrep(mdl.xvars,'_',' ')];
%             [~,inds] = sort([1:mdl.nx 1:mdl.nx]); legs = legs(inds); 
%             legend(legs);

            % plot summed coefficient timecourse
            subplot(2,1,2);
            plot(times,sum(scores_w_m,2),'k','linewidth',2); hold on;
            %
            xlim(lim.xlims.plot.LFP);
%             ylim([-0.4 0.4]); 
            set(gca,'YDir','reverse')
            plot(xlim,[0,0],'k-')
            xlabel 'time (s)'
            ylabel 'amplitude (uV)'

            %
            rs_saveExportFig(fig, figsdir, figname, true);
      
%         end

        %% PLOT FITTED EEG AVERAGE
        m = 1;
        mdl = models(m);

        % run estimated scores through mixing matrix to back-project to EEG channel-space
        scores_w_r = reshape(scores_w,ntrials*nsamps,ncomps);
        eeg_fitted_r = (scores_w_r / um');
        eeg_fitted = reshape(eeg_fitted_r,ntrials,nsamps,nchans);
        eeg_fitted_mean = squeeze(mean(eeg_fitted))';

%         % back-project to EEG channel-space (? equal to the version computed with single-trial scores)
%         eeg_fitted_mean2 = (scores_w_m / um')';

        % plot fitted EEG average
        figname = [ sub '_LME_LFP_EEG_PNT_EEG-FITTED_Cz_' s.conds{cond} ];
        fig = figure('name',figname, 'NumberTitle','off');
        inds = strcmp({ents_all.type},s.conds{cond});
        plotc = find(strcmp({data.EEG(1).chanlocs.labels},'CZ'));
        eeg_cond_mean = mean(deeg(:,:,inds),3);
        yyaxis 'left';  plot(times,  eeg_cond_mean(plotc,:)); ylabel 'EEG (uV)'; 
        set(gca,'YDir','reverse');
        xlim(lim.xlims.plot.EEG)
        yyaxis 'right'; plot(times,eeg_fitted_mean(plotc,:)); ylabel 'EEG fitted (A.U.)';
        xlim(lim.xlims.plot.EEG)

        rs_saveExportFig(fig, figsdir, figname);
%         close(fig)

        %% plot topographies

        figname = [ sub '_LME_LFP_EEG_PNT_TOPOMAX_' s.conds{cond} ]; 
        fig = figure('name',figname, 'NumberTitle','off');
        %
        vals = max(eeg_fitted_mean,[],2);
        ctemp = [];
        ctemp.clims = [-3 3]/10;
%         ctemp.clims = 'absmax'
        ctemp.colormap  = colormap2;
        ctemp.lay       = data.EEG(1).chanlocs;
        ctemp.gridscale = s.gridscale; % 700
        ctemp.numcontours = 0;
        rickoplot_MK( vals, ctemp );
        clim(ctemp.clims) % because the actual clims topoplot uses are slightly higher than requested ...
        colorbar

        rs_saveExportFig(fig, figsdir, figname);
%         close(fig)

        %% PLOT COMPONENT ESTIMATES AGAINST SUPRAMODALITY SCORE
        % ! wait instead of overall supramodality just compute selectivity for this particular condition!

        m = 1;
        mdl = models(m);

        figname = [ sub '_LME_LFP_EEG_PNT_SCAT-EST-SUP_' s.conds{cond} ]; 
        fig = figure('name',figname, 'NumberTitle','off');
        %
        est_abs = abs(mdl.est);
        
        % stats
%         [r,p] = corr(comp_sup,est_abs,'type','pearson');
        [r,p] = corr(comp_sup,est_abs,'type','spearman');

        % plot scatterplot
        scatter(comp_sup, est_abs); lsline; title(sprintf('r = %.2f  p = %.2f',r,p))

        % ? if doing stats here, need to consider that the supramodality metric range 
%            is from 1 to inf, while est_abs range is from 0 to inf


%         rs_saveExportFig(fig, figsdir, figname);
%         close(fig)


% %         % !! experimental version using only the predictors included in the best model
%         est_abs_INC  = abs(mdl.Coefficients.Estimate);
%         comp_sup_INC = comp_sup(selpreds);
%         comp_sup_EXC = comp_sup(setdiff(1:ncomps,selpreds));
% 
%         figure; scatter(comp_sup_INC, est_abs_INC); lsline
%         hold on; scatter(comp_sup_EXC, zeros(size(comp_sup_EXC))); lsline


   
        %%
        close all 
        

        %% CONDITION LOOP
    end 

%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
