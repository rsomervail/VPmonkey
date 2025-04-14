% 
% 
%       plot averages of LFP ONLY using LME estimates and confidence intervals 
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

% pool = gcp('nocreate');
% if isempty(pool)
%     pool = parpool(6);
% end

%% SETTINGS
s = [];
s.savePath =  [ getRoot '/VPmonkey/paper/results/lw' ];  % mkdir(s.savePath)

% plot limits 
lim = VPmonkey_fetchLimits;

% channels to plot
s.chans2plot.LFP = {'LFP'};

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

ftype = 'LFP';

DSF = 2; % downsample factor, 2 = 512, 4 = 256

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    savestr = []; % tracking which modules are used to reject trials

    %% LOAD DATA

    % by electrode (LFP, MUA etc)
    cfg = struct; cfg.exportlw = false;
    cfg.sub = sub;
    %

    % LFP electrode rejection on the basis of unusually high SD value
    cfg.ar_lfp_SD.thresh = 1.5;

%     cfg.ar_lfp.thresh_abs = [150,10/100; 300,0]; % default (reported in draft as of 20/03/2025)
% %     cfg.ar_lfp.thresh_abs = [150,5/100; 300,0]; % stricter version
%     cfg.ar_lfp.maxbad = 0.25; % default, rejects only one electrode (in SubT) (reported in draft as of 20/03/2025)
% %     cfg.ar_lfp.maxbad = 0.10; % stricter version that removes 4 electrodes total
    %
    cfg.filetypes = {ftype};
    cfg.LFP_include = {'LP_30'};
%     cfg.LFP_exclude = {'LP_30'};
    cfg.byElectrode = true;
    cfg.average   = false;
    cfg.mergesesh = false;
    %
    cfg.zscore = {};  % ? z-scoring may obscure important differences in VP presence!
%     cfg.zscore = {'LFP'};
%     cfg.zscore_win = [-0.2 0.6]; 
%     cfg.zscore_cond  = true;
    %
    [data,tbl_depths] = VPmonkey_mergeSesh(cfg);

    % downsample for speed
    data.(ftype) = pop_rs_downsample(data.(ftype),DSF);

    % extract time-window of interest
    data.(ftype) = pop_select(data.(ftype), 'time', lim.xlims.plot.(ftype));
    nelec = length(data.(ftype));

    %% LOOP THROUGH CONDITIONS
    for cond = 1:length(s.conds)
        close all % otherwise too many figures

        %% EXTRACT DATA CORRESPONDING TO CONDITION
       
        % get data from this condition for each electrode
        clear dcond ntrials_per_elec
        for e = 1:nelec
            ents_temp = data.LFP(e).event;
            inds = find(strcmp({ents_temp.type},s.conds{cond}));
            if ~isempty(inds)
                dcond(e) = pop_select(data.LFP(e), 'trial', inds );
            end
            ntrials_per_elec(e) = length(inds);
        end
    
        % extract all trials
        d = squeeze(cat(3,dcond.data));
        ents = cat(2,dcond.event);
        if length(ents) ~= size(d,2)
            error 'mismatch between events and ntrials - need to exclude repeat events'
        end
        [nsamps,ntrials] = size(d);
        d = d';
        times = dcond(1).times/1000;

        %% GET INFO FOR BUILDING PREDICTOR MATRICES

        % get per-trial session list
        seshlist = unique(tbl_depths.sesh);
        nsesh = length(seshlist);
        sesh = [];
        for e = 1:nelec
            sesh = [sesh; repmat( tbl_depths.sesh(e), ntrials_per_elec(e), 1) ];
        end

        % get per-trial depth values
        depthmet = 'Depth'; % ? can later try to see if there is any benefit to RelTOA or RelDura
        depth = [];
        for e = 1:nelec
            depth = [depth; ...
                repmat( tbl_depths.(depthmet)(e), ntrials_per_elec(e), 1) ];
        end
        depthz = (depth-mean(depth))/std(depth,1);

        % get per-trial binned depth values
        edges = [ min(depthz) quantile(depthz,2)  max(depthz) ];
        [counts,edges,bins] = histcounts(depthz,edges);
        depthBin_dummy = dummyvar(bins);
        depthBin_vars = arrayfun(@(x) ['depthBin_' num2str(x)], 1:length(counts), 'uniformoutput',false);

        % get per-trial XY coordinate values
        coordX = [];
        coordY = [];
        for e = 1:nelec
            coordX = [coordX; ...
                repmat( tbl_depths.x(e), ntrials_per_elec(e), 1) ];
            coordY = [coordY; ...
                repmat( tbl_depths.y(e), ntrials_per_elec(e), 1) ];
        end

%         % get depths relative to TOA   ? without missing values I can't compare this to other models
%         depthRelTOA = [];
%         for e = 1:nelec
%             depthRelTOA = [depthRelTOA; ...
%                 repmat( tbl_depths.DepthRelTOA(e), ntrials_per_elec(e), 1) ];
%         end
%         depthRelTOAz = (depthRelTOA-mean(depthRelTOA))/std(depthRelTOA,1);

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

        % random effects
        G = categorical(sesh);
        Z = ones(ntrials,1);

        %% PREPARE SET OF MODELS

        % build LME predictor matrices
        models = [];
        %
        models(end+1).name     = 'INT';
        models(end).formula    = 'y ~ 1 + (1 | s)'; %
        models(end).xvars      = {'Intercept'};
        models(end).x          = ones(ntrials,1);
        %
%         models(end+1).name    = 'INT_DEPTH';
%         models(end).formula = 'y ~ 1 + depth + (1 | s)'; %
%         models(end).x     = [ ones(ntrials,1)  depthz ]; 
%         models(end).xvars = {'Intercept','DepthZ'};
        %
        models(end+1).name    = 'INT_X';
        models(end).formula = 'y ~ 1 + Y + (1 | s)'; %
        models(end).x     = [ ones(ntrials,1) coordX  ]; 
        models(end).xvars = {'Intercept','X'};
        %
        models(end+1).name    = 'INT_Y';
        models(end).formula = 'y ~ 1 + Y + (1 | s)'; %
        models(end).x     = [ ones(ntrials,1)  coordY ]; 
        models(end).xvars = {'Intercept','Y'};
        %
        models(end+1).name    = 'INT_X_Y';
        models(end).formula = 'y ~ 1 + X + Y + (1 | s)'; %
        models(end).x     = [ ones(ntrials,1) coordX coordY ]; 
        models(end).xvars = {'Intercept','X','Y'};
        %
%         models(end+1).name    = 'INT_DEPTH_X';
%         models(end).formula = 'y ~ 1 + Y + (1 | s)'; %
%         models(end).x     = [ ones(ntrials,1) depthz coordX  ]; 
%         models(end).xvars = {'Intercept','DepthZ','X'};
%         %
%         models(end+1).name    = 'INT_DEPTH_Y';
%         models(end).formula = 'y ~ 1 + Y + (1 | s)'; %
%         models(end).x     = [ ones(ntrials,1) depthz  coordY ]; 
%         models(end).xvars = {'Intercept','DepthZ','Y'};
%         %
%         models(end+1).name    = 'INT_DEPTH_X_Y';
%         models(end).formula = 'y ~ 1 + X + Y + (1 | s)'; %
%         models(end).x     = [ ones(ntrials,1) depthz coordX coordY ]; 
%         models(end).xvars = {'Intercept','DepthZ','X','Y'};
        % 
        models(end+1).name    = 'SITE';
        models(end).formula = 'y ~ 1 + site + (1 | s)'; %
        models(end).x     = sitedummy; 
        models(end).xvars = sites';
%         %
%         models(end+1).name    = 'SITE_DEPTH';
%         models(end).formula = 'y ~ 1 + site + depth + (1 | s)'; %
%         models(end).x     = [ sitedummy depthz ]; 
%         models(end).xvars = [ sites', {'DepthZ'}];     
        %
        models(end+1).name    = 'DEPTHBIN';
        models(end).formula = 'y ~ depthBin_1 + depthBin_2 + ... + (1 | s)'; %
        models(end).x     = depthBin_dummy; 
        models(end).xvars = depthBin_vars;

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
                amp = double(d(:,t));
    
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
        resfile = [resdir filesep sub '_results_LME_LFPonly_' s.conds{cond} '.mat'];
        save(resfile,'models')

        %% (OPTIONAL) LOAD MODEL RESULTS
        if ~exist('res','var')
            load(resfile)
        end

        %% PLOT - MODEL BIC FITS (ALL)
        figname = [ sub '_LME_' ftype '_BIC_all_' s.conds{cond} ];

        % get BIC values
        models2plot = true(size(models));
        bic_all = cat(1,models(models2plot).bic);

        % compute relative BIC values for easier visualisation
%         temp2plot = bic_all - mean(bic_all); 
        temp2plot = bic_all - min(bic_all); % ? this is more interpretable than relative to mean

        % plot
        cols = distinguishable_colors(sum(models2plot));
        fig = figure('name',figname); 
        subplot(2,1,1);
        plot(times, mean(d),'LineWidth',2);
        subplot(2,1,2);
        for k = 1:sum(models2plot)
            plot(times, temp2plot(k,:), 'LineWidth',2, 'Color', cols(k,:)); hold on;
        end
        legend(strrep({models(models2plot).name},'_','-'));
        %
        rs_saveExportFig(fig, figsdir, figname, true);


        %% PLOT - MODEL BIC FITS (ONLY COORDS)
%         figname = [ sub '_LME_' ftype '_BIC_onlyCoords_' s.conds{cond} ];
% 
%         % get BIC values
%         models2plot = contains({models.name},'X') | contains({models.name},'Y');
%         bic_all = cat(1,models(models2plot).bic);
% 
%         % compute relative BIC values for easier visualisation
% %         temp2plot = bic_all; % sanity check 
% %         temp2plot = bic_all - mean(bic_all); 
%         temp2plot = bic_all - min(bic_all); % ? this is more interpretable than relative to mean
% 
%         % plot
%         cols = distinguishable_colors(sum(models2plot));
%         fig = figure('name',figname); 
%         subplot(2,1,1);
%         plot(times, mean(d),'LineWidth',2);
%         subplot(2,1,2);
%         for k = 1:sum(models2plot)
%             plot(times, temp2plot(k,:), 'LineWidth',2, 'Color', cols(k,:)); hold on;
%         end
%         legend(strrep({models(models2plot).name},'_','-'));
%         %
%         rs_saveExportFig(fig, figsdir, figname, true);

        %% PLOT - MODEL BIC FITS (NO COORDS)
%         figname = [ sub '_LME_' ftype '_BIC_noCoords_' s.conds{cond} ];
% 
%         % get BIC values
%         models2plot = ~(contains({models.name},'X') | contains({models.name},'Y'));
%         bic_all = cat(1,models(models2plot).bic);
% 
%         % compute relative BIC values for easier visualisation
% %         temp2plot = bic_all; % sanity check 
% %         temp2plot = bic_all - mean(bic_all); 
%         temp2plot = bic_all - min(bic_all); % ? this is more interpretable than relative to mean
% 
%         % plot
%         cols = distinguishable_colors(sum(models2plot));
%         fig = figure('name',figname); 
%         subplot(2,1,1);
%         plot(times, mean(d),'LineWidth',2);
%         subplot(2,1,2);
%         for k = 1:sum(models2plot)
%             plot(times, temp2plot(k,:), 'LineWidth',2, 'Color', cols(k,:)); hold on;
%         end
%         legend(strrep({models(models2plot).name},'_','-'));
%         %
%         rs_saveExportFig(fig, figsdir, figname, true);
        
        %% PLOT - MODEL BIC FITS (SHORTLIST OF TOP N MODELS WITHIN ERP TIMEWINDOW)
%         % ? ADDING TOGETHER RANKS LIKE THIS IS PRETTY WEIRD - not using for now
%         % ? TRYING to just add the subtracted versions although this is also weird
%         figname = [ sub '_LME_' ftype '_BIC_topN_' s.conds{cond} ];
% 
%         % define time-window and number of top models to plot
%         N = 4;
%         twin = findnearest(times,0):findnearest(times,0.3);
% 
%         % get BIC values
%         bic_all = cat(1,models.bic);
% 
% %         % (A) get models ranks within time-window   %!! DOESN'T ORDER PROPERLY
% %         [~,indy] = sort(bic_all);
% %         [~,indy] = sort(sum(indy(:,twin),2));
% %         models2plot(:) = false; models2plot(indy(1:N)) = true; 
% %         bic_all = cat(1,models(models2plot).bic);
% 
%         % compute relative BIC values for easier visualisation
% %         bic_all = bic_all - mean(bic_all); 
%         bic_all = bic_all - min(bic_all); % ? this is more interpretable than relative to mean
% 
%         % (B) compute summed relative BIC within window and use to rank models
%         [~,indy] = sort(sum(bic_all(:,twin),2));
%         indy = indy(1:N);
%         bic_all = bic_all(indy,:);
% 
%         % plot
%         cols = distinguishable_colors(sum(models2plot));
%         fig = figure('name',figname); 
%         subplot(2,1,1);
%         plot(times, mean(d),'LineWidth',2);
%         subplot(2,1,2);
%         for k = 1:N
%             plot(times, bic_all(k,:), 'LineWidth',2, 'Color', cols(k,:)); hold on;
%         end
%         legend(strrep({models(indy).name},'_','-'));
%         %
%         rs_saveExportFig(fig, figsdir, figname, true);

        %% LOOP THROUGH MODELS AND PLOT ALL COEFFICIENTS
        for m = 1:nmodels
            mdl = models(m);
            %
            figname = [ sub '_LME_' ftype '_COEF_' s.conds{cond} '_mdl-' mdl.name  ];

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
            if contains(mdl.name,'SITE')
                plot(times,mean(est2plot),'k','LineWidth',2)
            end
            xlim(lim.xlims.plot.(ftype));
            ylim(lim.ylims.plot.(ftype));
            set(gca,'YDir','reverse')
            plot(xlim,[0,0],'k-')
            xlabel 'time (s)'
            ylabel 'coefficient estimate (uV)'

            % legend
            legs = [repmat({''},1,mdl.nx), strrep(mdl.xvars,'_',' ')];
            [~,inds] = sort([1:mdl.nx 1:mdl.nx]);
            legs = legs(inds);
            if contains(mdl.name,'SITE')
                legs = [legs {''},{'AVG-SITE'}];
            end
            legend(legs);
            
            %
            rs_saveExportFig(fig, figsdir, figname, false);
      
        end

        %% PLOT - DEPTH - ELECTRODE AVGS COLOUR CODED BY DEPTH
        mdl = models(strcmp({models.name},'DEPTH'));
        if ~isempty(mdl)
            %
            figname = [ sub '_LME_' ftype '_DEPTHRANGE_' s.conds{cond} '_mdl-' mdl.name  ];

            % compute a range of depths
            ndivs = 100;
            depth2plot  =  linspace( min(depth), max(depth), ndivs )';
            depth2plotZ = (depth2plot-mean(depth))/std(depth,1)';

            % compute predicted mean responses in range
            int = mdl.est(1,:); % intercept
            dep = mdl.est(2,:); % depth estimates
            temp2plot = int + depth2plotZ * dep; % basically just the fixed effects part of the model formula
%             temp2plot = depth2plotZ * dep; % sanity check - just the depth coefficient
%             temp2plot = (1:ndeps)' * dep; % sanity check - linear increases

            % plot
            fig = figure('name',figname);
            cols = jet(ndivs); % blue means shallower, red means deeper
            for k = 1:ndivs
                plot(times, temp2plot(k,:), 'Color',cols(k,:), 'LineWidth',2); hold on;
            end
            % plot mean depth (i.e. intercept) as a thick black line
            plot(times, int, 'LineWidth',4,'Color','k'); 
            %
            xlim(lim.xlims.plot.(ftype));
            ylim(lim.ylims.plot.(ftype));
            set(gca,'YDir','reverse')
            plot(xlim,[0,0],'k-')
            xlabel 'time (s)'
            ylabel 'amplitude (uV)'

            %
            rs_saveExportFig(fig, figsdir, figname, false);

        end

        %% PLOT - X-Y - ELECTRODE AVGS COLOUR CODED BY X-Y
        mdl = models(strcmp({models.name},'INT_X_Y'));
        if ~isempty(mdl)
            %
            figname = [ sub '_LME_' ftype '_XYRANGE_' s.conds{cond} '_mdl-' mdl.name  ];

            % compute a range of coords
            ndivs = 100;
            x2plot  =  linspace( min(coordX), max(coordX), ndivs )';
            y2plot  =  linspace( min(coordY), max(coordY), ndivs )';

            % compute predicted mean responses in range
            int = mdl.est(1,:); % intercept
            x = mdl.est(2,:); % X estimates
            y = mdl.est(3,:); % Y estimates

            % average over X and Y respectively while varying the other
            temp2plotX = int + x2plot * x + mean(coordY) * y; 
            temp2plotY = int + y2plot * y + mean(coordX) * x; 

            % plot
            fig = figure('name',figname);
            cols = jet(ndivs); % blue means shallower, red means deeper
            for sbp = 1:2
                subplot(1,2,sbp)
                if sbp == 1
                    temp2plot = temp2plotX;
                    axtitle = 'varying X averaged over Y';
                else
                    temp2plot = temp2plotY;
                    axtitle = 'varying Y averaged over X';
                end

                %
                for k = 1:ndivs
                    plot(times, temp2plot(k,:), 'Color',cols(k,:), 'LineWidth',2); hold on;
                end
                % plot mean value (i.e. intercept) as a thick black line
                plot(times, int, 'LineWidth',4,'Color','k'); 
                %
                xlim(lim.xlims.plot.(ftype));
                ylim(lim.ylims.plot.(ftype));
                set(gca,'YDir','reverse')
                plot(xlim,[0,0],'k-')
                xlabel 'time (s)'
                ylabel 'amplitude (uV)'
                title(axtitle)
            end
            %
            rs_saveExportFig(fig, figsdir, figname, false);

        end

        %% CONDITION LOOP
    end 


%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
