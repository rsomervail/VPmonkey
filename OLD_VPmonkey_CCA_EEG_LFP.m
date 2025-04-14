%
%    run CCA between EEG and LFP, concatenating electrodes across all sessions as a single LFP channel
% 
%       
%       ? currently topography is ugly and scores seem to mostly reflect big high amplitude stuff 
%         in some middle trials/sessions (trial-wise mode)
%               - should add a thing to the spreadsheet & merge function for filtering bad LFP seshes
%       
%       ? MAYBE IT WOULD BE BETTER TO SWITCH THIS TO SOME KIND OF LINEAR REGRESSION WHICH TRIES TO MAXIMISE
%         REPLICATION OF THE DATA, RATHER THAN JUST MAXIMISING CORRELATION? SOME WAY TO DEAL WITH OVERFITTING TOO
%       
% 
%       OLD NOTES:
% 
%       ! should probs account for oversampling at middle depths .. how??
% 
%       ! BIGGEST issue seems to be that the topography seems to just reflect the intracranial recording hole
%           ? i.e. because there is a hole?
%           ? test this by doing a separate control CCA using only baseline (ideally totally outside of the VP epoch but that's more work)
% 
%       ! implement Sparse-CCA, find hyperparameters using cross-validation (it has a function for that) 
%           & use the cross-validation to check how regular the results are
% 
%       ? probs normalisation is not needed because CCA is looking for correlations which are scale invariant but worth doing anyway
%         for both LFP (done?) and EEG - note that this means I have to reverse z-scoring for any back-projected components
%           - normalisation should probs be done as described by Cohen GED tutorial: mean center each chan seperately, then divide by pooled SD (to preserve relationships across channels)
% 
%       ? check if PCA reduction is appropriate
%
%       ! also compare the strength of the CCA correlation with depth, both across whole thing and pbp
% 
%       ! can plot separate condition contribution timecourses + supramodality/selectivity timecourse
%           ! compute seperately for LFPs and test consistency of this
% 
%       ? computing CCA for the 3 conditions together avoids overfitting to a particular condition - but could arguably produce
%         the opposite bias? not sure, but SCCA should reduce overfitting anyhow
% 
% 
%       ! try adding big uncorrelated noise to data across a random set of channels and see if the CCA will
%           still recover the VP (tests to see if correlated comp should necessarily explain more of the avg variance)
%               ?? might need to make sure noise is phase-locked then?
%               ? the way I'm adding random noise now probably increases rank too
%               ? with super big random noise it was heavily reflected in the CCA...
% 
%       % ! subtract CC from data and export to lw for comparison ? but how to handle multiple LFPs/CCs in one dataset?
%           % ? dream result here would be removing basically everything from avg ERP except lemniscal components
% 
%       % ! maybe should compare the CC to permuted dataset CCA, either trials or all timepoints, or channels?
%           these all have different levels of severity and control for different things
% 
% 
%       % ? is CCA biased towards larger stuff in the EEG? (not because that stuff matches LFP better but just anyway)
%           ? adding noise seemed to imply this
%           ? could maybe test by normalising each timepoint 
% 
% 
%       ? warning: the subspace projection is not guaranteed to be correct for non-orthogonal components
%           so here I have the PCA (orthogonal) and CCA (not orthogonal, but only 1 comp so...?)
% 
%%
clc
clearvars
% close all

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

homedir = ['/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_' subfold];
cd(homedir);

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

cmap = rs_prepareTopo; % need to do this before loading eeglab

try 
    eeglab
catch
    addlab
end

%% SETTINGS
s = [];
s.savePath =  ['/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_' subfold]; % mkdir(s.savePath)

s.conds = {'VIS', 'SOM', 'AUD'}; nconds = length(s.conds);

% s.epwin = [ -0.2  0.7 ]; % window of interest
s.epwin = [ 0  0.5 ]; % ? doesn't allow comparing to baseline

s.lowpass = 30; % Hz, LFP, to match the EEG

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/neuro/iannettilab/Projects/VP monkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

s.filt_elec = '';
% s.filt_elec = 'Depth > MED'; % excluding median electrode here so that num electrodes is equal
% s.filt_elec = 'Depth < MED';
% s.filt_elec = 'DepthRelTOA > MED';
% s.filt_elec = 'DepthRelTOA < MED';

% SELECT ANIMAL
sub = 'SubM';
% sub = 'SubT';

%% import intracranial information table
tbl_depths = VPmonkey_import_electrodeDepths;

%% load data 

% by session (EEG)
cfg = struct;
cfg.sub = sub;
    % trial-rejection criteria
    ar = struct;
    ar.method = 'time';
    ar.metric = 'median';
    ar.thresh = 2; % <<<<<<<<<<<<<<<<<<<< AR SETTINGS ARE ESPECIALLY IMPORTANT FOR CCA ANALYSIS
    ar.timeprop = 0.15; % ?                 ! THERE ARE SOME REALLY BIG VALUES IN THE LFP CONCATED DATA
    ar.chanprop = 0.15; % ?
    cfg.autoar = ar;
cfg.filetypes = {'EEG','LFP'};
cfg.zscore    = {'EEG','LFP'};
cfg.zscore_cond = true;
cfg.byElectrode = true;
cfg.average   = true; % <<<<<<<<< AVERAGE TRIALS WITHIN SESSION/CONDITION OR NOT? BIG QUESTION
cfg.mergesesh = true;
cfg.filt_elec = s.filt_elec;
cfg.LFP_include = {'LP_30'};
cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
cfg.EEG_exclude = {'noICA'};
cfg.filt_preproc = 'allGood';  % these refer to EEG preproc

[data, tbl_depths] = VPmonkey_mergeSesh(cfg);
%
nfiles = size(tbl_depths,1);
tbl_depths.ntrials = data.EEG.etc.merge.ntrials;
%
trialinds = [];
for k = 1:nfiles
    trialinds = [trialinds; k * ones(tbl_depths.ntrials(k),1) ];
end

% CROP ALL DATA TO A COMMON TIMEWINDOW (ASSUMES SRATE IS EQUIVALENT)
data.EEG = pop_select(data.EEG, 'time', [-0.2 0.6] ); % ? MAYBE INCLUDES TOO MUCH BASELINE?
data.LFP = pop_select(data.LFP, 'time', [-0.2 0.6] );

%% LOOP THROUGH CONDITIONS
ents = data.EEG.event;
for cond = 1:length(s.conds)

    %% EXTRACT DATA AND CONCATENATE
    indy = strcmp({ents.type},s.conds{cond});
    tbl_depths2 = tbl_depths(trialinds(indy),:);

    deeg  = double(data.EEG.data(:,:,indy));
    dlfp  = double(data.LFP.data(:,:,indy));
    
    [nchans,nsamps,ntrials]  = size(deeg);
    
    deegr = reshape(deeg, nchans, nsamps * ntrials)';
    dlfpr = reshape(dlfp,      1, nsamps * ntrials)';

    % center
    deegr = deegr - mean(deegr);
    dlfpr = dlfpr - mean(dlfpr);
    
    %% CANONICAL CORRELATION ANALYSIS
%     close all
%     
%     % RUN CCA
%     disp 'running CCA ...'
%     tic
% 
%     reg_eeg = 0.01; % ? L1 regularisation is for when you think few variables have non-zero weights ...
%     reg_lfp = 0.01; %   ? not sure if this is what I need to reduce overfitting to certain observations
%     [coef_eeg,coef_lfp,rvals,scores_eeg_r,scores_lfp_r] = scanoncorr(deegr,dlfpr,reg_eeg,reg_lfp);
% 
% %     [coef_eeg,coef_lfp,rvals,scores_eeg,scores_lfp] = canoncorr(deegr,dlfpr);
% 
%     toc
%     disp '... CCA finished'
%     % ! make EEGLAB implementation of scanoncorr which is nicer formatted
%     
% 
%     % ANALYSIS OF CCA RESULTS
%     close all
% 
%     times = data.EEG.times/1000;
% 
% %     % plot EEG coefficients as a topoplot
% %     figure; rs_topoplot(coef_eeg,chanlocs);
% %     %! export topoplot w/ nicer Brewer colourmap
% 
%     % scale EEG data timepoints by corresponding scores and plot mean topography)
%     scores_eeg_r_norm = scores_eeg_r/max(abs(scores_eeg_r));
%     scores_eeg_r_norm .* deegr;
%     deegr_scaled = scores_eeg_r_norm .* deegr;
%     figure; rs_topoplot(mean(deegr_scaled)',chanlocs); %! export topoplot w/ nicer Brewer colourmap
% 
%     % explore scores per trial/avg
%     scores_eeg = reshape(scores_eeg_r,nsamps,ntrials);
%     scores_lfp = reshape(scores_lfp_r,nsamps,ntrials);
%     figure; 
%     subplot(3,1,1); plot(times, scores_eeg); hold on; plot(times,mean(scores_eeg,2),'k','LineWidth',3); title 'scores - EEG';
%     subplot(3,1,2); plot(times, scores_lfp); hold on; plot(times,mean(scores_lfp,2),'k','LineWidth',3); title 'scores - LFP';
%     subplot(3,1,3); yyaxis left; plot(times, mean(scores_eeg,2)); yyaxis right; plot(times,mean(scores_lfp,2)); title 'scores - BOTH';
% 
%     % ! look at scores to determine things about lfp locations/depths/etc
%    
    %% LINEAR MODEL - PREDICTING EEG FROM LFP
    
%     % timepoints as observations
%     mdls = cell(1,nchans);
%     topo = nan(nchans,1);
%     fitted_r = nan(nsamps*ntrials,nchans);
%     for c = 1:nchans
%         mdls{c}     = fitlm( dlfpr, deegr(:,c) );
%         topo(c,1)   = mdls{c}.Coefficients.Estimate(2);
%         pvals(c,1)  = mdls{c}.Coefficients.pValue(2);
%         fitted_r(:,c) = mdls{c}.Fitted;
%     end
%     fitted_r = fitted_r';
%     fitted = reshape(fitted_r,nchans,nsamps,ntrials);
%     figure; rs_topoplot(topo,chanlocs); %! export topoplot w/ nicer Brewer colourmap

    % ? maybe go back to running this for every session separately?

    % timepoints as separate models
    mdls = cell(nchans,nsamps);
%     fitted = nan(nchans,nsamps,ntrials);
    dlfp2 = squeeze(dlfp)';
    deeg2 = permute(deeg,[3 1 2]);
    parfor c = 1:nchans
        mdls_temp = cell(1,nsamps);
        deeg2_temp = squeeze(deeg2(:,c,:));
        for t = 1:nsamps
            mdls_temp{1,t} = fitlm( dlfp2(:,t), deeg2_temp(:,t) );
%             fitted(c,t,:) = mdls{c,t}.Fitted;
        end
        mdls(c,:) = mdls_temp;
        %
        c
    end

    % ? ADD LINEAR OR NON-LINEAR (POWERLAW?) EFFECT OF TRIAL NUMBER

    % ? CONSIDER USING MULTIVARIATE REGRESSION - right now I effectively have an intercept per timepoint in the epoch...

    %% extract info from models
    fitted     = nan(nchans,nsamps,ntrials);
    bic        = nan(nchans,nsamps);
    rsq        = nan(nchans,nsamps);
    est_int_mean = nan(nchans,nsamps);
    est_lfp_mean = nan(nchans,nsamps);
    est_int_se = nan(nchans,nsamps);
    est_lfp_se = nan(nchans,nsamps);
    est_int_ci  = nan(nchans,nsamps,2);
    est_lfp_ci  = nan(nchans,nsamps,2);
    est_cov    = nan(nchans,nsamps); 
    for c = 1:nchans
        for t = 1:nsamps
            fitted(c,t,:) = mdls{c,t}.Fitted';
            bic(c,t)      = mdls{c,t}.ModelCriterion.BIC;
            rsq(c,t)      = mdls{c,t}.Rsquared.Ordinary;
            est_int_mean(c,t) = mdls{c,t}.Coefficients.Estimate(1);
            est_lfp_mean(c,t) = mdls{c,t}.Coefficients.Estimate(2);
            est_int_se(c,t) = mdls{c,t}.Coefficients.SE(1); % standard error of intercept coefficient
            est_lfp_se(c,t) = mdls{c,t}.Coefficients.SE(2); % standard error of LFP predictor coefficient estimate 
            est_cov(c,t) = mdls{c,t}.CoefficientCovariance(1,2); % covariance between intercept (mean value) and weight of LFP predictor coefficient
            ci = mdls{c,t}.coefCI;
            est_int_ci(c,t,:) = ci(1,:); % confidence intervals (? 95%) for intercept 
            est_lfp_ci(c,t,:) = ci(2,:); % confidence intervals (? 95%) for LFP predictor 
        end
    end

    times = data.EEG.times/1000;
    t1 = findnearest(times,0);
    t2 = findnearest(times,0.3);

%     residuals = deeg - fitted;
%     res_sum = squeeze(sum(abs( residuals(:,t1:t2,:) ), 1:2));

    %% plot
    close all

    % single-channel model summary plot
    figure('name','model parameters over time'); 
    plotchan = find(strcmpi({data.EEG.chanlocs.labels},'cz'));
    xlims = [min(times),max(times)];
    nrows = 5; 

    subplot(nrows,1,1);
    plot(times, mean(deeg(plotchan,:,:),3)); hold on;
    plot(times, mean(dlfp,3));
    plot(xlims,[0,0],'k--')
    xlabel 'time (s)'; ylabel 'EEG amplitude (uV)'; 
    title 'EEG vs LFP'; 
    legend({'EEG','LFP'})
    xlim(xlims)

    subplot(nrows,1,2);
%     plot(times,est_int_mean(plotchan,:)); hold on
    bounds = squeeze( abs( est_int_ci(plotchan,:,:) - est_int_mean(plotchan,:) ));  % CIs are actually ~ 2xSE - wonder if I should just use that?
    boundedline(times',est_int_mean(plotchan,:), bounds,'cmap',[0,0,1]); hold on
%     plot(times,est_lfp_mean(plotchan,:));
    bounds = squeeze( abs( est_lfp_ci(plotchan,:,:) - est_lfp_mean(plotchan,:) ));  % CIs are actually ~ 2xSE - wonder if I should just use that?
    boundedline(times,est_lfp_mean(plotchan,:), bounds,'cmap',[1,0,0]); hold on
    plot(xlims,[0,0],'k--')
    xlim(xlims)
    legend({'INT_S_E','INT_M_E_A_N','LFP_S_E','LFP_M_E_A_N'})
    title 'estimate mean'

    subplot(nrows,1,3);
    plot(times,est_int_se(plotchan,:)); hold on
    plot(times,est_lfp_se(plotchan,:));
    plot(xlims,[0,0],'k--')
    xlim(xlims)
    legend({'INT','LFP'})
    title 'estimate SE'

    % BIC
    subplot(nrows,1,4);
    plot(times, bic(plotchan,:))
    xlim(xlims)
    title 'BIC'
    % lowest BIC is best

    % rsquared
    subplot(nrows,1,5);
    % unscaled
    yyaxis left
    plot(times, rsq(plotchan,:)); hold on; 
    % scaled by variance
    yyaxis right
    rsq_norm_var = rsq .* var(deeg,[],3); % ? should this be population variance? since we are interested in the sample itself?
    plot(times, rsq_norm_var(plotchan,:)); 
%     % scaled by mean
%     rsq_norm_mean = rsq .* abs(mean(deeg,3));
%     plot(times, rsq_norm_mean(plotchan,:)); 
    %
    xlim(xlims)
    legend({'unscaled','norm-var'})
    title 'R-squared'
    
    %% topoplot model values   
    %? could plot a nice long series of these at regular intervals, matching colour scales
    %   - with rows for parameters/variables
    lat = 0.09375; % peak of BIC (highest = worst)
    t = findnearest(times,lat);
    figure('name',sprintf('topoplots - t = %.2f',lat)); 
    ncols = 10;

    subplot(1,ncols,1); 
    topo = est_int_mean(:,t);
    rs_topoplot(topo,chanlocs); %! export topoplot w/ nicer Brewer colourmap    
    title 'est - INT - mean'

    subplot(1,ncols,2); 
%     topo = est_int_se(:,t); 
    topo = est_int_se(:,t) ./ var(deeg(:,t,:),[],3);
    rs_topoplot(topo,chanlocs); %! export topoplot w/ nicer Brewer colourmap    
%     title 'est - INT - SE'
    title 'est - INT - SE norm-var'

    subplot(1,ncols,3); 
    topo = est_lfp_mean(:,t);
    rs_topoplot(topo,chanlocs); %! export topoplot w/ nicer Brewer colourmap    
    title 'est - LFP - mean'

    subplot(1,ncols,4); 
%     topo = est_lfp_se(:,t);
    topo = est_lfp_se(:,t) ./ var(deeg(:,t,:),[],3);
    rs_topoplot(topo,chanlocs); %! export topoplot w/ nicer Brewer colourmap    
%     title 'est - LFP - SE'
    title 'est - LFP - SE norm-var'
    
    subplot(1,ncols,5); 
    topo = bic(:,t);
    rs_topoplot(topo,chanlocs); %! export topoplot w/ nicer Brewer colourmap    
    title 'BIC'  % ? IS IT REASONABLE TO COMPARE BIC BETWEEN MODELS PREDICTING DIFFERENT DATA?
    
    subplot(1,ncols,6);
    topo = rsq(:,t);
    rs_topoplot(topo,chanlocs); %! export topoplot w/ nicer Brewer colourmap    
    title 'rsq';
%     subplot(1,ncols,7);
%     topo = rsq_norm_mean(:,t);
%     rs_topoplot(topo,chanlocs); %! export topoplot w/ nicer Brewer colourmap    
%     title 'rsq-norm-mean';
    subplot(1,ncols,7);
    topo = rsq_norm_var(:,t); % ? in the case of topoplots, normalising per channel may obscure topographic differences
    rs_topoplot(topo,chanlocs); %! export topoplot w/ nicer Brewer colourmap    
    title 'rsq-norm-var';


    % ! should REALLY check the SE and such across channels too, e.g. plotting them across scalp at relevant timepoints


% ?? maybe would be good to compare these metrics vs an intercept only model? (feels like nullism tho...)
    %? would R2 increase? BIC?


%% END LOOP THROUGH CONDS
end


%% END
% cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


