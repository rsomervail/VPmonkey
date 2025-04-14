% 
% 
%       plot averages of EYE only using LME estimates and confidence intervals 
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

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

s.conds = {'AUD','SOM','VIS'}; nconds = length(s.conds);

subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};

DSF = 2; % downsample factor, 2 = 512, 4 = 256

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};

    %% LOAD DATA

%     % trial-rejection criteria
%     ar = struct;
%     ar.method = 'time';
%     ar.metric = 'median';
%     ar.thresh = 3;
%     ar.timeprop = 0.1; % ? too strict?
%     ar.chanprop = 0.1; % ? too strict?

    cfg = struct; cfg.exportlw = false;
    cfg.sub = sub;
%     cfg.ar = ar;
    cfg.filetypes = {'EYE'};
%     cfg.filt_eye = true; % only get good EYE sessions
    cfg.filt_eye = false; % ? KEEPING ALL EYE SESSIONS 
    cfg.EYE_include = {'LP_5'}; % only include low-pass filtered pupil data
%     cfg.EYE_exclude = {'LP_5'}; % exclude low-pass filtered pupil data for GED filtering, do LP here instead if necessary
    cfg.byElectrode = false;
    cfg.average   = false;
    cfg.mergesesh = false;

%     % z-score
%     cfg.zscore = {'EYE'};  
%     cfg.zscore_win.EYE = lim.xlims.plot.EYE;
%     cfg.zscore_cond  = false; % ? here not z-scoring per condition, to allow qualtitative comparisons between conditions
    
    data = VPmonkey_mergeSesh(cfg);
    nsesh = length(data.EYE);

    % downsample for speed
    if DSF > 1
        data.EYE = pop_rs_downsample(data.EYE,DSF);
    end

    % extract time-window of interest
    data.EYE = pop_select(data.EYE, 'time', lim.xlims.plot.EYE);
    

    %% LOOP THROUGH CONDITIONS
    for cond = 1:length(s.conds)

        %% EXTRACT DATA CORRESPONDING TO CONDITION
       
        % get data for this condition from each session 
        clear dcond ntrials_per_sesh
        for sh = 1:nsesh
            ents_temp = data.EYE(sh).event;
            inds = find(strcmp({ents_temp.type},s.conds{cond}));
            if ~isempty(inds)
                dcond(sh) = pop_select(data.EYE(sh), 'trial', inds );
            end
            ntrials_per_sesh(sh) = length(inds);
        end
    
        % extract all trials
        deye = squeeze(cat(3,dcond.data))';
        ents = cat(2,dcond.event);
        [ntrials, nsamps] = size(deye);
        if length(ents) ~= ntrials
            error 'mismatch between events and ntrials - need to exclude repeat events'
        end
        times = dcond(1).times/1000;

        %% FIX EYE ARTIFACTS WITH TEMPORAL PCA/GED FILTER (? UNUSED IN THE END DUE TO REMAINING ARTIFACT BEING "SMEARED OUT")
% 
%         % OPTIONAL - plot data before filtering 
%         figure; 
%         subplot(2,2,1); plot(times,deye); tempylims = ylim; % ? IMPORTANT - look at remaining artifacts and think how to model
%         subplot(2,2,2); plot(times,mean(deye));
%         
%         % PCA version
%         temp = double(deye);
%         [coef, scores, evals, ~, expvar, mu] = pca(temp); 
% %         c = 40; figure; plot(coef(:,c)); % plot a coefficient
% %         ncomps = 3; % 3-4 keeps most of the variance
%         ncomps = find( cumsum(expvar)>90, 1, 'first'); % ? this might vary with downsample rate...
%         sum(expvar(1:ncomps))
%         deye = (mu + scores(:,1:ncomps) * coef(:,1:ncomps)');
% 
% %         % GED version - BROKEN CURRENTLY
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
%         % plot PCA-cleaned data
%         subplot(2,2,3); plot(times,deye); ylim(tempylims) % ? IMPORTANT - look at remaining artifacts and think how to model
%         subplot(2,2,4); plot(times,mean(deye))

        %% GET INFO FOR BUILDING PREDICTOR MATRICES

        % get per-trial session list 
        sesh = [];
        for sh = 1:nsesh
            sesh = [sesh; repmat( sh, ntrials_per_sesh(sh), 1) ];
        end 
      
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
       
        % get number of variables for each model
        nmodels = length(models); for m = 1:nmodels, models(m).nx = length(models(m).xvars); end

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
        resfile = [resdir filesep sub '_results_LME_EYEonly_' s.conds{cond} '.mat'];
        save(resfile,'models')

        %% (OPTIONAL) LOAD MODEL RESULTS
        if ~exist('res','var')
            load(resfile)
        end

        %% PLOT - MODEL BIC FITS
%         figname = [ sub '_LME_' 'EYE' '_BIC_' s.conds{cond} ];
% 
%         % compute relative BIC values for easier visualisation
%         bic_all = cat(1,models.bic);
% %         temp2plot = bic_all - mean(bic_all); 
%         temp2plot = bic_all - min(bic_all); % ? this is more interpretable than relative to mean
% 
%         % plot
%         fig = figure('name',figname); 
%         subplot(2,1,1);
%         plot(times, mean(deye),'LineWidth',2);
%         subplot(2,1,2);
%         plot(times, temp2plot, 'LineWidth',2);
%         legend(strrep({models.name},'_','-'));
%         %
%         rs_saveExportFig(fig, figsdir, figname, true);

        %% LOOP THROUGH MODELS AND PLOT ALL COEFFICIENTS
        for m = 1:nmodels
            mdl = models(m);
            %
            figname = [ sub '_LME_EYE_COEF_' s.conds{cond}  ];
%             figname = [ sub '_LME_EYE_COEF_' s.conds{cond} '_mdl-' mdl.name  ];

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
            xlim(lim.xlims.plot.EYE);
            ylim(lim.ylims.plot.EYE); 
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
        close all

        %% CONDITION LOOP
    end 


%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
