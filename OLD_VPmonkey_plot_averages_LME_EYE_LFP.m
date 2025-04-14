% 
% 
%       plot averages of EYE using LME estimates and confidence intervals 
%           - this version also includes LFP as a predictor to demonstrate relationship between pupil response & VP
%           
%           - this old version replicated all the EYE trials once per electrode... 
%               replaced with per-session version using LFP pooled electrode average in new version
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

% plot limits  !! MAKE THIS INTO A GLOBAL FUNCTION!!
% xlims
s.xlims.plot.EEG  = [-0.2 0.6];  
s.xlims.plot.LFP  = [-0.2 0.6];  
s.xlims.plot.MUA  = [-0.2 0.6];  
s.xlims.plot.EYE  = [-0.5   3]; 
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

split_range  = [0, 0.4]; % range to compute ERP amplitude within

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

s.conds = {'AUD','SOM','VIS'}; nconds = length(s.conds);

% subs = {'SubM','SubT'};
% subs = {'SubM'};
subs = {'SubT'};

DSF = 2; % downsample factor, 2 = 512, 4 = 256

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
%     ar.thresh = 3; % ? lower if not much difference
%     ar.timeprop = 0.1; % 
%     ar.chanprop = 0; % one channel in EYE data
%     ar.filetypes = {'EYE'}; % ? LFP best handled by absolute here because I don't want to regress to the mean

    % ar_LFP  % ? should eventually replace this with a more general ar_abs setting
    ar_lfp.thresh_abs = [150,10/100; 300,0]; % default (reported in draft as of 20/03/2025)
%     ar_lfp.thresh_abs = [150,5/100; 300,0]; % stricter version
    ar_lfp.maxbad = 0.25; % default, rejects only one electrode (in SubT) (reported in draft as of 20/03/2025)
%     ar_lfp.maxbad = 0.10; % stricter version that removes 4 electrodes total

    % ar_EYE
    ar_eye.thresh_sign = [-0.05,3/100];
    ar_eye.timewin = [-0.5 3];
    ar_eye.maxbad  = 0.25;

    % by electrode (LFP, MUA etc)
    cfg = struct; cfg.exportlw = false;
    cfg.sub = sub;
    if exist('ar','var'), cfg.ar = ar; end
    if exist('ar_lfp','var'), cfg.ar_lfp = ar_lfp; end
    if exist('ar_eye','var'), cfg.ar_eye = ar_eye; end
    cfg.filetypes = {'LFP','EYE'};
    cfg.LFP_include = {'LP_30'};
%     cfg.LFP_exclude = {'LP_30'};
%     cfg.filt_eye = true; % ONLY KEEP GOOD EYE sessions
    cfg.filt_eye = false; % ? KEEP ALL EYE SESSIONS - actually better SNR this way!
    cfg.EYE_include = {'LP_5'}; % only include low-pass filtered pupil data
    cfg.byElectrode = true;
    cfg.average   = false;
    cfg.mergesesh = false;

    % z-score
    cfg.zscore = {'LFP'};  % ? here z-scoring seems helpful for improving interpretation of coefficients (also because overall LFP amplitude varies with electrode location/depth)
%     cfg.zscore = {'EYE','LFP'}; % ? not doing EYE for consistency with LME_EYE
%     cfg.zscore_win.EYE = s.xlims.plot.EYE;
    cfg.zscore_win.LFP = s.xlims.plot.LFP;
    cfg.zscore_cond   = true;

    [data,tbl_depths] = VPmonkey_mergeSesh(cfg);

    % downsample for speed
    data.LFP = pop_rs_downsample(data.LFP,DSF/2);
    data.EYE = pop_rs_downsample(data.EYE,DSF);

    % extract time-window of interest
    data.LFP = pop_select(data.LFP, 'time', s.xlims.plot.LFP);
    data.EYE = pop_select(data.EYE, 'time', s.xlims.plot.EYE);
    nelec = length(data.LFP);

    %% LOOP THROUGH CONDITIONS
    for cond = 1:length(s.conds)

        %% EXTRACT DATA CORRESPONDING TO CONDITION
       
        %
        clear dcond
        for e = 1:nelec
            ents_temp = data.LFP(e).event;
            inds = find(strcmp({ents_temp.type},s.conds{cond}));
            if ~isempty(inds)
                dcond.LFP(e) = pop_select(data.LFP(e), 'trial', inds );
                dcond.EYE(e) = pop_select(data.EYE(e), 'trial', inds );
            end
            ntrials_per_elec(e) = length(inds);
        end

        % extract all trials
        deye = squeeze(cat(3,dcond.EYE.data))';
        dlfp = squeeze(cat(3,dcond.LFP.data))';
        ents = cat(2,dcond.LFP.event);
        if length(ents) ~= size(deye,1)
            error 'mismatch between events and ntrials - need to exclude repeat events'
        end
        [ntrials,nsamps] = size(deye);
        times = dcond.EYE(1).times/1000;

%         % OPTIONAL - plot all trials 
%         figure; plot(times,deye) % ? IMPORTANT - look at remaining artifacts and think how to model

        %% FIX EYE ARTIFACTS WITH TEMPORAL PCA/GED FILTER (? UNUSED IN THE END DUE TO REMAINING ARTIFACT BEING "SMEARED OUT")
        % ! flipped the deye matrix, needs to be edited to use this section

%         % PCA version   ? should really be done with non-replicated trials and before low-pass filter
%         temp = double(deye)';
%         [coef, scores, evals, ~, expvar, mu] = pca(temp); 
% %         c = 3; figure; plot(coef(:,c)); % plot a coefficient
%         ncomps = 3; % 3-4 keeps most of the variance
%         sum(expvar(1:ncomps))
%         deye = (mu + scores(:,1:ncomps) * coef(:,1:ncomps)');
% 
% %         % GED version - ? couldn't get it to work...
% %         temp = double(deye);
% %         mu = mean(temp,2);
% % %         temp_c = temp - mu; % ? is this centering done correctly??
% %         temp_c = temp; % ! sanity check with no centering
% %         [nvars,ntrials] = size(temp);
% %         ctemp = [];
% %         ctemp.dataR = temp_c;
% %         ctemp.dataS = repmat( mean(temp_c,2),1,ntrials);
% %         [evals, inverse, forward, scores] = rs_ged(temp_c, ctemp);
% 
% 
%         % plot PCA-cleaned trials
%         figure; plot(times, deye)
%         figure; plot(times, mean(deye))

        %% GET INFO FOR BUILDING PREDICTOR MATRICES

        % get per-trial session list
        seshlist = unique(tbl_depths.sesh);
        nsesh = length(seshlist);
        sesh = [];
        for e = 1:nelec
            sesh = [sesh; repmat( tbl_depths.sesh(e), ntrials_per_elec(e), 1) ];
        end

        % get per-trial depth values
        depthmet = 'Depth'; % ? can later try to see if there is any benefit to RelTOA or RelaDura
        depth = [];
        for e = 1:nelec
            depth = [depth; ...
                repmat( tbl_depths.(depthmet)(e), ntrials_per_elec(e), 1) ];
        end
        depthz = (depth-mean(depth))/std(depth,1);

        % get per-trial site values
        site = [];
        for e = 1:nelec
            site = [site; ...
                repmat( tbl_depths.dist_angle(e), ntrials_per_elec(e), 1) ];
        end
        [sites,~,sitenum] = unique(site);
        sites = strrep(strrep(sites,'-','m'),'.','d');
        sites = cellfun(@(x) ['site_' x], sites, 'UniformOutput',false);
        sitedummy = dummyvar(sitenum);

        % get per trial LFP amplitude
        times_lfp = dcond.LFP(1).times/1000;
        inds = findnearest(times_lfp,split_range(1)):findnearest(times_lfp,split_range(2));
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
        models(end).formula = 'eye ~ 1 + (1 | s)'; %
        % 
%         models(end+1).name  = 'INT_LFPrankz';
%         models(end).x     = [ ones(ntrials,1) lfp_rankz ]; 
%         models(end).xvars = {'Intercept','LFPrankz'};
%         models(end).formula = 'eye ~ 1 + lfp_rankz + (1 | s)'; %
        % 
%         models(end+1).name  = 'INT_LFPlog';
%         models(end).x     = [ ones(ntrials,1) log(lfp) ]; 
%         models(end).xvars = {'Intercept','LFPlog'};
%         models(end).formula = 'eye ~ 1 + log(lfp) + (1 | s)'; %
%         % 
%         models(end+1).name  = 'INT_LFPz';
%         models(end).x     = [ ones(ntrials,1) lfpz ]; 
%         models(end).xvars = {'Intercept','LFPz'};
%         models(end).formula = 'eye ~ 1 + lfpz + (1 | s)'; %
        % 
%         models(end+1).name  = 'LFPz';
%         models(end).x     = [ lfpz ]; 
%         models(end).xvars = {'LFPz'};
%         models(end).formula = 'eye ~ lfpz + (1 | s)'; %
%         % 
%         models(end+1).name    = 'INT_DEPTHxLFPz';
%         models(end).x     = [ ones(ntrials,1) depthz.*lfpz ];  
%         models(end).xvars = {'Intercept', [depthmet 'ZxLFPz']};
%         models(end).formula = 'eye ~ 1 + depth:lfpz + (1 | s)'; %        
%         %
%         models(end+1).name    = 'DEPTHxLFP_LFP';
%         models(end).x     = [ lfpz,  depthz.*lfpz ]; 
%         models(end).xvars = {'LFPz',[depthmet 'ZxLFPz']};
%         models(end).formula = 'eye ~ depth:lfpz + lfpz + (1 | s)'; %
%         % 
        models(end+1).name    = 'LFPz_bin';
        models(end).x     = [ lfpz_bin ]; 
        models(end).xvars = arrayfun(@(x) ['LFPz_bin_' num2str(x)] , 1:nbins, 'UniformOutput',false);
        models(end).formula = 'eye ~ lfpz_bin + (1 | s)'; %
        %
%         models(end+1).name    = 'DEPTHxLFPz_bin';
%         models(end).x     = [ depthz .* lfpz_bin ]; 
%         models(end).xvars = arrayfun(@(x) [depthmet 'ZxLFPz_bin_' num2str(x)] , 1:nbins, 'UniformOutput',false);
%         models(end).formula = 'eye ~ lfpz_bin + (1 | s)'; %
%         %
%         models(end+1).name    = 'LFPz_bin_DEPTHxLFPz_bin';
%         models(end).x     = [ lfpz_bin, depthz .* lfpz_bin ]; 
%         xvar_bin = arrayfun(@(x) ['LFPz_bin_' num2str(x)] , 1:nbins, 'UniformOutput',false);
%         xvar_binXdepth = arrayfun(@(x) [depthmet 'ZxLFPz_bin_' num2str(x)] , 1:nbins, 'UniformOutput',false);
%         models(end).xvars = [xvar_bin, xvar_binXdepth];
%         models(end).formula = 'eye ~ lfpz_bin + (1 | s)'; %

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
                amp = double(deye(:,t)); 
    
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
        resfile = [resdir filesep sub '_results_LME_EYE_LFP_' s.conds{cond} '.mat'];
        save(resfile,'models')

        %% (OPTIONAL) LOAD MODEL RESULTS
        if ~exist('models','var')
            load(resfile)
        end

        %% PLOT - MODEL BIC FITS
        figname = [ sub '_LME_EYE_LFP_BIC_' s.conds{cond} ];

        % compute relative BIC values for easier visualisation
        bic_all = cat(1,models.bic);
%         bic_all = bic_all - mean(bic_all); 
        bic_all = bic_all - min(bic_all); % ? this is more interpretable than relative to mean

        % plot
        fig = figure('name',figname); 
        subplot(2,1,1);
        plot(times, mean(deye),'LineWidth',2);
        subplot(2,1,2);
        plot(times, bic_all, 'LineWidth',2);
        legend(strrep({models.name},'_','-'));
        %
        rs_saveExportFig(fig, figsdir, figname, true);

        %% LOOP THROUGH MODELS AND PLOT ALL COEFFICIENTS
        for m = 1:nmodels
            mdl = models(m);
            %
            figname = [ sub '_LME_EYE_LFP_COEF_' s.conds{cond} '_mdl-' mdl.name  ];

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
            xlim(s.xlims.plot.EYE);
%             ylim(s.ylims.plot.EYE); % ? not sure what to set limits to yet
%             set(gca,'YDir','reverse')
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
%         close all 
        

        %% CONDITION LOOP
    end 

%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
