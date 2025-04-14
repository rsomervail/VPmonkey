% 
% 
%       plot averages of all data types across all sessions
%       plot t-test results 
% 
% 
%  
%       - Richard Somervail, 2023
%%
clc
clearvars
close all

homedir = '/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_main';
cd(homedir);

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

%% SETTINGS
s = [];
figsdir = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/paper/figures/raw';
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
s.ylims.plot.LFP  = [-4 4];   
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

% file_types = {'EYE','EEG','LFP','MUA'};
% file_types = {'EYE','EEG'};
% file_types = {'EYE','EEG','LFP'};
file_types = {'EYE'};
% file_types = {'EEG'};
% file_types = {'LFP'};


%% get topoplot stuff
addpath([  getRoot filesep 'MATLAB' filesep 'cbrewer' ])

[colormap2]=cbrewer('div', 'RdBu', 1000, 'linear'); % cubic , linear
colormap2 = colormap2(length(colormap2):-1:1, : );

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    savestr = []; % tracking which modules are used to reject trials

    %% LOAD DATA

    % trial-rejection criteria
    ar = struct;
    ar.method = 'time';
    ar.metric = 'median';
    ar.thresh = 3;
    ar.timeprop = 0.1; % ? too strict?
    ar.chanprop = 0.1; % ? too strict?

    % by session (EEG, EYE)
    cfg = struct;
    cfg.sub = sub;
    cfg.autoar = ar;
    cfg.byElectrode = false;
    cfg.average   = true;
    cfg.mergesesh = true;
    % EEG
    if any(ismember(file_types,{'EEG'}))
        cfg.filetypes = {'EEG'};
        cfg.filt_preproc = 'allGoodEEG'; % 'allGoodEEG', 'allPerfEEG', 'allSesh'
        cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
        cfg.EEG_exclude = {'noICA'};
        
        data_EEG = VPmonkey_mergeSesh(cfg);
        data.EEG = data_EEG.EEG
    end
    % EYE
    if any(ismember(file_types,{'EYE'}))
        cfg.filetypes = {'EYE'};
        cfg.filt_eye = true;
        cfg.EYE_include = {'LP_5'};
        %
        data_EYE = VPmonkey_mergeSesh(cfg);
        data.EYE = data_EYE.EYE;
    end

    % by electrode (LFP, MUA etc)
    if any(ismember(file_types,{'MUA','LFP'}))
        cfg = struct;
        cfg.sub = sub;
        cfg.autoar = ar;
        cfg.filetypes = file_types(ismember(file_types,{'MUA','LFP'}));
        cfg.LFP_include = {'LP_30'};
        cfg.byElectrode = true;
        cfg.average   = true;
        cfg.mergesesh = true;
%         cfg.zscore = [];
        cfg.zscore = cfg.filetypes;  % here I am using z-score because it's just a simple average and can't handle differences in overall amplitude
        cfg.zscore_win = [-0.2 0.6];
        data_elec = VPmonkey_mergeSesh(cfg);

        for k = 1:length(cfg.filetypes)
            data.(cfg.filetypes{k}) = data_elec.(cfg.filetypes{k});
        end
    end
    
    clearvars data_* 

    %% loop through file_types 
    for ft = 1:length(file_types)
        ftype = file_types{ft};

        %% get events
        dtemp = data.(ftype);
        t0s = findnearest(dtemp.times,0);
        t0s = t0s : dtemp.pnts : (dtemp.trials-1) * dtemp.pnts  + t0s;
        ents = dtemp.event;
        inds = ismember([ents.latency], t0s);
        ents = ents(inds);

        %% loop through conditions
        for cond = 1:length(s.conds)

            %% split by condition
            dcond = pop_select( data.(ftype), ...
                'trial', find(strcmp({ents.type},s.conds{cond}))  );

            %% select only the time-window to be plotted
            dcond = pop_select(dcond, 'time', s.xlims.plot.(ftype) );
            times = dcond.times/1000;
            xlims = s.xlims.plot.(ftype); % for actual setting of xlimits

            %% run cluster-permutation test
            ctemp = struct;
            ctemp.nperms = 1000;
    
%             ctemp.alpha1 = 0.05; % could lower this threshold if no clusters are formed, without affecting false positive rate
            ctemp.alpha1 = 0.1; % decided to lower to improve forming of clusters. can cite Maris & Oostenweld to justify
            ctemp.alpha2 = 0.05; % alpha-level for cluster significance - do not change
    
%             clocal.cluster_metric = 'maxT';  % fewer false positives
            ctemp.cluster_metric = 'meanT'; % fewer false negatives  ? is this compatible with lower alpha1 threshold?
    
            ctemp.minnumchan = 0; % ? for now allowing single-channel clusters and hoping their cluster-T will be low
    %                                  if there are many spurious single-channel clusters, then increase this
    
            if dcond.nbchan > 1
                ctemp.neighbours = 0.8;  % produces an average of ~9 channels
            end
            
            % RUN CLUSTER PERMUTATION ONE-SAMPLED T-TEST
            [dcond_tvals,dcond_pvals,dcond_thresh, dcond_clusterids] = pop_rs_clusterperm_ttest(dcond,ctemp);

%             % OPTIONAL PLOT
%             pop_letsplot([dcond_tvals, dcond_pvals dcond_thresh, dcond_clusterids])

            % EXPORT RESULTS TO LW
            cd([homedir filesep 'lw'])
            if strcmp(ftype,'EEG')
                templocs = chanlocs_lw;
            else
                templocs = [];
            end
            condstr = ['ep ' s.conds{cond} ' '];
            rs_convert_lab2lw_V1( dcond_tvals,      [ condstr dcond_tvals.setname],      templocs );
            rs_convert_lab2lw_V1( dcond_pvals,      [ condstr dcond_pvals.setname],      templocs );
            rs_convert_lab2lw_V1( dcond_thresh,     [ condstr dcond_thresh.setname],     templocs );
            rs_convert_lab2lw_V1( dcond_clusterids, [ condstr dcond_clusterids.setname], templocs );

            %% select only subset of channels to be plotted
            % data
            dcond = pop_select(dcond,'channel', ...
                find(ismember({dcond.chanlocs.labels},s.chans2plot.(ftype) )) );
            dcond.setname = ['ep ' s.conds{cond} ' ' dcond.setname];
            plotchans = {dcond.chanlocs.labels};

            % t-test results
            dcond_tvals = pop_select(dcond_tvals,'channel', ...
                find(ismember({dcond_tvals.chanlocs.labels},s.chans2plot.(ftype) )) );
            dcond_pvals = pop_select(dcond_pvals,'channel', ...
                find(ismember({dcond_pvals.chanlocs.labels},s.chans2plot.(ftype) )) );

            %% loop through channels to plot
            for c = 1:length(plotchans)

                %% get this channel
                dtemp = squeeze( dcond.data(c,:,:) );
                dtemp = permute(dtemp,[2 1]); % put epoch dimension first for later mean & t-test

                %% extract t-value and p-value timecourses for cluster permutation t-test for this channel
                tvals           = squeeze(dcond_tvals.data(c,:,1));
                tvals_thresh    = squeeze(dcond_tvals.data(c,:,2)); 
                pvals           = squeeze(dcond_pvals.data(c,:,1));
                pvals_thresh    = squeeze(dcond_pvals.data(c,:,2)); 

                %% take average of this channel
                davg = mean(dtemp, 'omitnan');
                davg_thresh = davg;
                davg_thresh(pvals_thresh==1) = 0;

                %% PLOT - average
                figname = [ sub '_avg_'  ftype '_' plotchans{c} '_' s.conds{cond} ]; 
                figs.(figname) = figure('name',figname, 'NumberTitle','off');

                % plot
                plot(times, davg, 'k', 'LineWidth',1 ); hold on;
                plot(times, davg_thresh, 'k', 'LineWidth',2.5 ); 
                xlim(xlims);
                ylim(s.ylims.plot.(ftype));
                plot(xlim,[0,0],'k')

                %% PLOT - average w/ bounded lines for Confidence Intervals
%                 figname = [ sub '_avgCI_'  ftype '_' plotchans{c} '_' s.conds{cond} ]; 
%                 figs.(figname) = figure('name',figname, 'NumberTitle','off');
%                 % plot
%                 bounds = double(CI - davg)'; % ?? because CI are absolute points either side of mean but function wants distance from davg
%                 [hl, hp] = boundedline( times, double(davg), abs(bounds), 'k'); hold on;
%                 xlim(xlims);
%                 ylim(s.ylims.plot.(ftype));
%                 plot(xlim,[0,0],'k')
              
                %% PLOT - t-values
                figname = [ sub '_tvals_'  ftype '_' plotchans{c} '_' s.conds{cond} ]; 
                figs.(figname) = figure('name',figname, 'NumberTitle','off');
                % plot
                plot(times, tvals, 'b' ); hold on;
%                 plot(times, tvals_thresh, 'k', 'LineWidth',2 );  % ? no plotting formed clusters atm
                xlim(xlims);
                ylim(s.ylims.tvals.(ftype));

                %% PLOT - p-values
%                 figname = [ sub '_pvals_'  ftype '_' plotchans{c} '_' s.conds{cond} ]; 
%                 figs.(figname) = figure('name',figname, 'NumberTitle','off');
%                 % plot
%                 plot(times, pvals, 'b' ); hold on;
% %                 plot(times, pvals_thresh, 'k', 'LineWidth',2 );  % ? no plotting formed clusters atm
%                 xlim(xlims);
%                 ylim([0,1]);

                %% PLOT - p-values (surprisal)
                figname = [ sub '_psurp_'  ftype '_' plotchans{c} '_' s.conds{cond} ]; 
                figs.(figname) = figure('name',figname, 'NumberTitle','off');
                % plot
                plot(times, -log2(pvals), 'b' ); hold on;
%                 plot(times, -log2(pvals_thresh), 'k', 'LineWidth',2 );  % ? no plotting formed clusters atm
                xlim(xlims);
                ylim([0,20]);
                plot(xlims,repmat(-log2(0.05),1,2),'k:')
                ylabel 'surprisal (-log2(p))'

            %% END loop through channels
            end

            %% PLOT - EEG topographies
            if strcmp(ftype,'EEG')

                % get copy of EEG data
                dcond = pop_select( data.(ftype), ...
                'trial', find(strcmp({ents.type},s.conds{cond}))  );
                dcond = pop_select(dcond, 'time', s.xlims.plot.(ftype) );
                times = dcond.times/1000;
                
                % get peaks to plot for this condition
                desired_peaklats = s.topopeaks.(sub).(s.conds{cond});
                npeaks = length(desired_peaklats);

                % find the closest real peaks on the average to this point
                peakchan = find(strcmp({dcond.chanlocs.labels},'CZ'));
                davg = mean(dcond.data,3); % make average of the peak channel
                [~,all_peaklocs] = findpeaks(abs(davg(peakchan,:))); 
                all_peaklats = times(all_peaklocs);
                final_lats = nan(1,npeaks);
                final_locs = nan(1,npeaks);
                for pk = 1:npeaks
                    final_lats(pk) = all_peaklats( findnearest( all_peaklats , desired_peaklats(pk) ) );
                    [~, final_locs(pk)] = min(abs( times - final_lats(pk) ));
                end

                % plot vertical lines on EEG average to indicate where the final topos are
                figure( figs.([ sub '_avg_'  ftype '_' plotchans{c} '_' s.conds{cond} ])  )
                for pk = 1:npeaks
                    plot([final_lats(pk),final_lats(pk)],ylim,'k')
                end

                % loop through identified peaks and plot topographies
                for pk = 1:npeaks
            
                    figname = [ sub '_topo_'  ftype '_' s.conds{cond} '_' num2str(pk) ]; 
                    figs.(figname) = figure('name',figname, 'NumberTitle','off');
                    % plot
                    vals = davg(:,final_locs(pk));
                    ctemp = [];
                    ctemp.clims     = s.clims.(sub).(s.conds{cond});
                    ctemp.colormap  = colormap2;
                    ctemp.lay       = dcond.chanlocs;
                    ctemp.gridscale = 150; % 700
                    ctemp.numcontours = 0;
                    rickoplot_MK( vals, ctemp );

                end
            end

            %% SAVE & EXPORT ALL CURRENT FIGURES
            cd(figsdir)
            fignames  = fields(figs);
            nfigs = length(fignames);
            for f = 1:nfigs
                rs_saveExportFig(figs.(fignames{f}), figsdir, fignames{f});
            end
            clear figs
            close all

        %% END loop through conditions
        end
    end
    %% END loop through file_types

%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
