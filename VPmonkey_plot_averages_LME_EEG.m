% 
% 
%       plot averages of EEG ONLY using LME estimates and confidence intervals 
% 
%       - Richard Somervail, 2023
%%
clc
clearvars
close all

homedir = '/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_main';
cd(homedir);

figsdir =  [ getRoot '/VPmonkey/paper/figures/raw' ]; 
resdir =  [ getRoot '/VPmonkey/paper/results' ]; 

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

% EEG topo peaks  ? this will find the closest peak on the average to the specified points
s.topopeaks.SubM.AUD = [30 95 160]/1000; 
s.topopeaks.SubM.SOM = [36 97 177] / 1000;
s.topopeaks.SubM.VIS = [22 43  67 98 120 210 290 ]/1000;
s.topopeaks.SubT.AUD = [33 79 150]/1000;
s.topopeaks.SubT.SOM = [36 120 170] / 1000;
s.topopeaks.SubT.VIS = [19 46  110  190 290 ]/1000;

% channels to plot
s.chans2plot.EEG = {'CZ'};

s.conds = {'AUD','SOM','VIS'}; nconds = length(s.conds);

subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};

file_types = {'EEG'};

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

    % get data by session 
    cfg = struct;
    cfg.sub = sub;
    cfg.byElectrode = false;
    cfg.average   = false;
    cfg.mergesesh = false;
    %
    cfg.filetypes = {'EEG'};
    cfg.filt_preproc = 'allSesh'; % 'allGood' or 'allSesh'  
    cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
    cfg.EEG_exclude = {'noICA'};
    %
    data = VPmonkey_mergeSesh(cfg);
    nsesh = length(data.EEG);
   
    % downsample for speed
    data.EEG = pop_rs_downsample(data.EEG,DSF);

    % extract time-window of interest
    data.EEG = pop_select(data.EEG, 'time', lim.xlims.plot.EEG);

    % get channel info
    chans = {data.EEG(1).chanlocs.labels};

    %% LOOP THROUGH CONDITIONS
    for cond = 1:length(s.conds)

        %% EXTRACT DATA CORRESPONDING TO CONDITION

        % get data for this condition from each session 
        clear dcond ntrials_per_sesh
        for sh = 1:nsesh
            ents_temp = data.EEG(sh).event;
            inds = find(strcmp({ents_temp.type},s.conds{cond}));
            if ~isempty(inds)
                dcond(sh) = pop_select(data.EEG(sh), 'trial', inds );
            end
            ntrials_per_sesh(sh) = length(inds);
        end

        % extract all trials
        d = squeeze(cat(3,dcond.data));
        ents = cat(2,dcond.event);
        [nchans,nsamps,ntrials] = size(d);
        if length(ents) ~= ntrials
            error 'mismatch between events and ntrials - need to exclude repeat events'
        end
        times = dcond(1).times/1000;

        % EXTRACT DATA FOR CHANNEL
        plotchan = find(strcmpi(chans,s.chans2plot.EEG{1}));
        dchan    = double( squeeze(d(plotchan,:,:))' );

        %% GET INFO FOR BUILDING PREDICTOR MATRICES

        % get per-trial session list 
        sesh = [];
        for sh = 1:nsesh
            sesh = [sesh; repmat( sh, ntrials_per_sesh(sh), 1) ];
        end 

        % fixed effects
        X = ones(ntrials,1); % intercept predictor for this condition

        % random effects
        G = categorical(sesh);
        Z = ones(ntrials,1);

        %% PREPARE MODEL OUTPUTS
        est = nan(1,nsamps);
        upp = nan(1,nsamps);

        %% LOOP THROUGH TIMEPOINTS
        fprintf('%s',repmat('.',1,nsamps))
        fprintf('\n\n')
        tin = tic;
        parfor t = 1:nsamps
%             t = 304; 
            amp = dchan(:,t);
   
            %% FIT MODEL

            % fit linear mixed-effects model
            mdl = fitlmematrix( X, amp, Z, G, ...
                'FitMethod','REML');

            % store outputs
            est(1,t)   = mdl.Coefficients.Estimate;
            upp(1,t)   = mdl.Coefficients.Upper;
        
            fprintf('\b|\n');
%             fprintf('|');
        %% TIME LOOP
        end
        fprintf('LME fitting complete in %.2f mins\n\n',toc(tin)/60)

        %% COMPUTE MARGIN OF ERROR (RELATIVE CONFIDENCE INTERVAL)
        moe = upp - est;

        % STORE OUTPUTS
        models = struct; % for now only one model is fit for EEG
        models.est = est;
        models.moe = upp - est; % COMPUTE MARGIN OF ERROR (RELATIVE CONFIDENCE INTERVAL)

        %% SAVE MODEL RESULTS
        resfile = [resdir filesep sub '_results_LME_EEGonly_' s.conds{cond} '.mat'];
        save(resfile,'models')

        %% (OPTIONAL) LOAD MODEL RESULTS
        if ~exist('res','var')
            load(resfile)
        end

        %% PLOT AT CZ
        
        % plot mean coefficient estimate w/ confidence intervals
        figname = [ sub '_LME_EEG_COEF_' s.chans2plot.EEG{1} '_' s.conds{cond} ];
        fig = figure('name',figname); 
        boundedline(times, est, moe );
        hold on;
        xlim(lim.xlims.plot.EEG);
        ylim(lim.ylims.plot.EEG);
        plot(xlim,[0,0],'k-')
        xlabel 'time (s)'
        ylabel 'coefficient estimate (uV)'
        set(gca,'YDir','reverse')
        %
        rs_saveExportFig(fig, figsdir, figname);
            
        %% PLOT TOPOGRAPHIES 

        % get peaks to plot for this condition
        desired_peaklats = s.topopeaks.(sub).(s.conds{cond});
        npeaks = length(desired_peaklats);

        % find the closest real peaks on the average to this point
        [~,all_peaklocs] = findpeaks(abs(est)); 
        all_peaklats = times(all_peaklocs);
        final_lats = nan(1,npeaks);
        final_locs = nan(1,npeaks);
        for pk = 1:npeaks
            final_lats(pk) = all_peaklats( findnearest( all_peaklats , desired_peaklats(pk) ) );
            [~, final_locs(pk)] = min(abs( times - final_lats(pk) ));
        end

        % run LME for all channels for each chosen timepoint 
        est_chans = nan(nchans,npeaks);
        for pk = 1:npeaks
            for c = 1:nchans
                amp = double( squeeze(d(c,final_locs(pk),:)) );
                mdl = fitlmematrix( X, amp, Z, G, 'FitMethod','REML');
                est_chans(c,pk) = mdl.Coefficients.Estimate;
            end
        end

        % plot vertical lines on EEG average to indicate where the final topos are
        figname = [ sub '_LME_EEG_topoPEAKS_' s.chans2plot.EEG{1} '_' s.conds{cond} ];
        fig = figure;
        plot(times,est); hold on;
        for pk = 1:npeaks
            plot([final_lats(pk),final_lats(pk)],ylim,'k')
        end
        rs_saveExportFig(fig, figsdir, figname);
        close(fig)

        % loop through identified peaks and plot topographies
        cd(figsdir)
        for pk = 1:npeaks
            figname = [ sub '_LME_EEG_topo_' s.conds{cond} '_' num2str(pk) ]; 
            fig = figure('name',figname, 'NumberTitle','off');

            % plot
            vals = est_chans(:,pk);
            ctemp = [];
%             ctemp.clims     = lim.clims.(sub).(s.conds{cond});
                clims = ceil(max(abs(vals))); % like absmax but nicely rounded
                clims = [-clims clims];
            ctemp.clims     = clims;
            ctemp.colormap  = colormap2;
            ctemp.lay       = data.EEG(1).chanlocs;
            ctemp.gridscale = s.gridscale; % 700
            ctemp.numcontours = 0;
            rickoplot_MK( vals, ctemp );
            clim(clims) % because the actual clims topoplot uses are slightly higher than requested ...
            colorbar

            rs_saveExportFig(fig, figsdir, figname);
            close(fig)

        end
        
        %% CONDITION LOOP
    end 
    close all

%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
