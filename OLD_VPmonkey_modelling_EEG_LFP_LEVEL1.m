%
% 
% 
%    Run linear models between EEG and LFP, concatenating electrodes across all sessions as a single LFP channel
% 
% 
%           !! ADD COMPONENT METRICS, E.G. SUPRAMODALITY, WIDESPREADNESS, DISTANCE OF MASS FROM ELECTRODE CHAMBER ETC
%               - could do all this in second level script but then would need to reload EEG data or just use scores
%
% 
%           !!  BUG: MISMATCH BETWEEN FILE 22 TRIALS AND EVENTS FOR COND 3 AND MAYBE THE OTHERS?
%               - different number of EEG trials in different repeats of the same data... 
%                   ** because of trial rejection? !! **
%               - testing by deactivating AR for now, hopefully handled by robust fitting and COMP analysis done on averages
%  
%           !! once pipeline is finished, re-run with robust fitting!!!
% 
%            !! then can plot model fit Rsqa against component properties, electrode depths etc
%                as second order model, but using SE of underlying predictor estimates?
%               ! make sure to include possible confound of component overall variance when assessing supramodality & widespreadness
%                   ? is it actually a confound or just inevitable?
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

s.epwin = [ -0.2  0.6 ]; % window of interest
% s.epwin = [ 0  0.5 ]; 

s.lowpass = 30; % Hz, LFP, to match the EEG

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/neuro/iannettilab/Projects/VP monkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

[colormap2]=cbrewer('div', 'RdBu', 1000, 'linear'); % cubic linear
colormap2 = colormap2(length(colormap2):-1:1, : );

% SELECT ANIMAL
sub = 'SubM';
% sub = 'SubT';

%% import intracranial information table
tbl_depths = VPmonkey_import_electrodeDepths;

%% load data 

% by session (EEG)
cfg = struct;
cfg.sub = sub;
%     % trial-rejection criteria   ? switched off to try to prevent trial-mismatch bug on file 22 (SubM)
%     ar = struct;
%     ar.method = 'time';
%     ar.metric = 'median';
%     ar.thresh = 2; 
%     ar.timeprop = 0.15; % ?                
%     ar.chanprop = 0.15; % ?
%     cfg.autoar = ar;
cfg.filetypes = {'EEG','LFP'};
cfg.zscore    = {'EEG','LFP'};
cfg.zscore_cond = true;
cfg.byElectrode = 1; % treat each LFP channel as one channel over time & repeat EEG
% cfg.byElectrode = 2; % treat each LFP channel as a separate channel, which is nan in trials outside the session it belongs too
% cfg.average   = true; % <<<<<<<<< AVERAGE TRIALS WITHIN SESSION/CONDITION OR NOT? BIG QUESTION
% cfg.mergesesh = true;
cfg.average   = false; % <<<<<<<<< AVERAGE TRIALS WITHIN SESSION/CONDITION OR NOT? BIG QUESTION
cfg.mergesesh = false;
cfg.filt_elec = '';
cfg.LFP_include = {'LP_30'};
cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
cfg.EEG_exclude = {'noICA'};
cfg.filt_preproc = 'allGood';  % these refer to EEG preproc

[data, tbl_depths] = VPmonkey_mergeSesh(cfg);
%
% ntrials_perSesh = data.EEG.etc.merge.ntrials;
% nfiles  = length(ntrials_perSesh);
% 
% %
% trialinds = [];
% for k = 1:nfiles
%     trialinds = [trialinds; k * ones(ntrials_perSesh(k),1) ];
% end

% CROP ALL DATA TO A COMMON TIMEWINDOW (ASSUMES SRATE IS EQUIVALENT)
data.EEG = pop_select(data.EEG, 'time', s.epwin ); 
data.LFP = pop_select(data.LFP, 'time', s.epwin );

% DOWNSAMPLE
for f = 1:length(data.EEG)
    data.EEG(f) = pop_rs_downsample(data.EEG(f), 2);
    data.LFP(f) = pop_rs_downsample(data.LFP(f), 2);
end
% get times
times = data.EEG(1).times/1000;

%% COMPONENT ANALYSIS

% merge trials across sessions
seshlist = unique(tbl_depths.sesh);
nsesh = length(seshlist);
deeg_avg = cell(nsesh,1);
ntrials_eeg_persesh = nan(nsesh,1);
for sh = 1:nsesh
    indy = find(strcmp(tbl_depths.sesh,seshlist{sh}),1,'first');
    deeg_avg{sh} = mean(data.EEG(indy).data,3);
    ntrials_eeg_persesh(sh) = size(data.EEG(indy).data, 3);
end
deeg_avg = cat(3,deeg_avg{:});

% prepare data for component analysis
[nchans,nsamps,~] = size(deeg_avg);
deeg_avg_r = reshape(deeg_avg,nchans,nsamps*nsesh)';
mu_eeg = mean(deeg_avg_r);
deeg_avg_r  = deeg_avg_r - mu_eeg;

% COMPONENT ANALYSIS
% PCA
[pca_coef, pca_scores_r, pca_latent] = pca(deeg_avg_r);

% account for rank deficiency
comps2keep = pca_latent > 1e-7;
ncomps2keep = sum(comps2keep);
pca_coef     = pca_coef(:,comps2keep);
pca_scores_r = pca_scores_r(:,comps2keep)';
pca_scores = reshape(pca_scores_r, ncomps2keep, nsamps, nsesh); 

% run ICA
[ica_weights,ica_sphere,ica_compvars,~,~,~,ica_act_r] = runica(deeg_avg_r','pca',ncomps2keep);
ica_winv = pinv(ica_weights*ica_sphere); % pseudo inverse to get forward (?) weights for topo-plotting
ica_ncomps = length(ica_compvars);
ica_act = reshape(ica_act_r, ica_ncomps, nsamps, nsesh);

% run GED 
%!

% select component type to use:
%     comptype = 'PCA';
comptype = 'ICA';
%     comptype = 'GED';
%

% get components to use for the model
switch comptype
    case 'PCA'
        coef2use     = pca_coef;  
        coef2use_inv = pca_coef';
        scores2use   = pca_scores;
    case 'ICA'
        coef2use     = ica_winv;
        coef2use_inv = ica_weights*ica_sphere;
        scores2use   = ica_act;
    case 'GED'
        %!!!!!
end
ncomps = size(coef2use,2);

% loop through components and plot topographies 
topo = coef2use;
figure('name','component topographies');
for k = 1:ncomps
    subplot(5,5,k);
    rs_topoplot(topo(:,k), chanlocs, 'colormap',colormap2); 
    title(num2str(k))
%         c = 5; clim([-c,c])
end

% apply spatial filter to single trials and compute scores
dcomp = cell(nsesh,1);
for sh = 1:nsesh
    indy   = find(strcmp(tbl_depths.sesh,seshlist{sh}),1,'first');
    temp   = data.EEG(indy).data;
    temp_r = reshape(temp, nchans,nsamps*size(temp,3));
    temp   = reshape(coef2use_inv*temp_r,ncomps,nsamps,size(temp,3));
    dcomp{sh} = temp;
end

% compute metric of "widespread-ness" for each component topography 
wideness = nan(ncomps,1);
for c = 1:ncomps
    topo = coef2use(:,c);
    cov = topo * topo'; % compute covariance matrix
    % normalise to population pearson coefficient
    for row = 1:nchans 
        for col = 1:nchans
            cov(row,col) = cov(row,col) ./ (sqrt(cov(row,row))*sqrt(cov(col,col)));
        end
    end
    % sum of absolute correlations, excluding self-correlation, divided by perfect correlation of all electrodes
    wideness(c) = (sum(abs(cov),'all')-nchans) / nchans^2; 
end

% !! COMPUTE COMPONENTS SUPRAMODALITY AND WIDESPREADNESS SCORES
% ! USE compvar FUNCTION ON PER CONDITION SET OF TRIALS WITHIN 0 TO 0.5 WINDOW TO QUANTIFY SUPRAMODALITY

% ! then need to store this information with results structure


% save memory by clearing EEG data
data = struct('LFP',data.LFP);

%% LOOP THROUGH CONDS
nfiles = length(data.LFP);
est    = nan(nconds,nfiles,ncomps,nsamps);
SE     = nan(nconds,nfiles,ncomps,nsamps);
rsq    = nan(nconds,nfiles,ncomps,nsamps); 
for cond = 1:nconds
%         cond = 3; warning 'SKIPPING CONDITIONS 1 & 2'  % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
    %% LOOP THROUGH ELECTRODES (FILES)
    for f = 1:nfiles
%         f = 22; warning 'FOCUSING ON FILE 22 BUG'  % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< DEBUG <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        % GET SESSION DATA
        seshnum = find(strcmp(tbl_depths.sesh{f}, seshlist));
        deeg = dcomp{seshnum};
        dlfp = squeeze(data.LFP(f).data);
    
        % GET CONDITION
        ents = {data.LFP(f).event.type};
        indy = strcmp(ents,s.conds{cond});
        
        %% IF THIS CONDITION HAS TRIALS IN THIS FILE, THEN CONTINUE, OTHERWISE SKIP FILE
        if any(indy)
            deeg  = deeg(:,:,indy);
            dlfp2 = dlfp(:,indy);

            %% LOOP THROUGH COMPONENTS
            for c = 1:ncomps
                deeg2 = squeeze(deeg(c,:,:));
    
                %% LOOP THROUGH TIMEPOINTS
                est_temp = nan(nsamps,1);
                SE_temp  = nan(nsamps,1);
                rsq_temp = nan(nsamps,1);
                parfor t = 1:nsamps
    
                    % prepare data for model fitting
                    X = dlfp2(t,:);
                    Y = deeg2(t,:);
                    
                    % standardise data  % <<< does this mean the information about mean EEG inside LFP is removed?
                    X = (X - mean(X)) ./ std(X);
                    Y = (Y - mean(Y)) ./ std(Y);
        
                    % fit model
                    m = fitlm(X,Y,'intercept',false);  % >>>> USE ROBUST FITTING BUT INCREASE ITERATION LIMIT
%                     m = fitlm(X,Y,'intercept',false, 'RobustOpts','on' );  % ? using default robust fitting method
                    est_temp(t) = m.Coefficients.Estimate;
                    SE_temp(t)  = m.Coefficients.SE;
                    rsq_temp(t) = m.Rsquared.Ordinary;
    
                %% END LOOP THROUGH TIMEPOINTS
                end
                
                % store values from these models
                est(cond,f,c,:) = est_temp;
                SE(cond,f,c,:)  = SE_temp;
                rsq(cond,f,c,:) = rsq_temp;
    
            %% END LOOP THROUGH COMPONENTS
            end

        %% END CHECK FOR TRIALS WITHIN CONDITION
        end

    %% END FILE LOOP
    fprintf('file %03d/%03d\n',f,nfiles)
    end

%% END LOOP THROUGH CONDS
end

%% SAVE RESULTS
% rename outputs
LFP_est = est;
LFP_SE  = SE;
%
varlist = {'LFP_est','LFP_SE','rsq','tbl_depths','coef2use','wideness','times','s.conds'}; 
r = struct; for k = 1:length(varlist), results.(varlist{k}) = eval(varlist{k}); end
save([getRoot '/VPmonkey/paper/results/'  sub '_results_LM_LFP_EEG_LEVEL1.mat'], 'r', '-v7.3') 


%% COMPUTE SECOND ORDER MODEL DO AS SECOND SCRIPT WHICH CALLS R SCRIPT EVENTUALLY

% !! SPLIT THIS OFF INTO A SECOND LEVEL SCRIPT & RENAME THIS ONE ACCORDINGLY

% ! try various predictors and predict both rsquared and coefficient and compare their BIC values
% ? plot a range of models, not just single best BIC per timepoint/condition/whatever

% ! have separate models for looking at condition effects, component variance/supramodality, electrode depths & locations?

% ? really would be ideal to factor in uncertainty (SE) in coefficient estimates directly by multi-level modelling somehow
%       ? maybe should do this bit in R  (? can call a script from here)


%% TEMP PLOTTING
cond = 3;
rsqa2plot = squeeze(rsq(cond,:,:,:));

rsqa2plot   = rsqa2plot(1:f-1,:,:);
tbl_depths2 = tbl_depths(1:f-1,:);

rsqa2plot = atanh(rsqa2plot); % Fisher transform

t = findnearest(times,0.1);
% t = findnearest(times,0.2);

rsqa2plot = rsqa2plot(:,:,t);

% loop through components and plot mean rsqa per component against electrode depth
figure;
metric = 'DepthRelTOA';
for c = 1:ncomps
    subplot(5,5,c);
    scatter(tbl_depths2.(metric), rsqa2plot(:,c));
    axis square
end

% plot component Rsqa as a boxplot  ? violin would be nicer eventually
figure;
boxplot(rsqa2plot);

%% OLD 
% %% LOOP THROUGH CONDITIONS
% ents = data.EEG.event;
% % for cond = 1:length(s.conds)
% for cond = 3:length(s.conds) % <<<<<<<<<<<<<<<<<<<<<<<, DEBUG
% 
%     %% EXTRACT DATA AND CONCATENATE
%     indy = strcmp({ents.type},s.conds{cond});
%     trialinds2  = trialinds(indy);
%     tbl_depths2 = tbl_depths(trialinds2,:);
% 
%     deeg  = double(data.EEG.data(:,:,indy));
%     dlfp  = double(data.LFP.data(:,:,indy));
%     
%     [nchans,nsamps,ntrials]  = size(deeg);
%     
%     % reshape & center
%     deegr  = reshape(deeg, nchans, nsamps * ntrials)';
%     dlfpr  = reshape(dlfp,       1, nsamps * ntrials)';
%     mu_eeg = mean(deegr);
%     deegr  = deegr - mu_eeg;
%     mu_lfp = mean(dlfpr,'omitnan');
%     dlfpr  = dlfpr - mu_lfp;  % ? not sure about centering, does it push everything up or down?   
% 
%     %% COMPONENT ANALYSIS
% 
%     % PCA
%     [pca_coef, pca_scores_r, pca_latent] = pca(deegr);
%     
%     % account for rank deficiency
%     comps2keep = pca_latent > 1e-7;
%     ncomps2keep = sum(comps2keep);
%     pca_coef     = pca_coef(:,comps2keep);
%     pca_scores_r = pca_scores_r(:,comps2keep)';
%     pca_scores = reshape(pca_scores_r, ncomps2keep, nsamps, ntrials); 
% 
%     % run ICA
%     [ica_weights,ica_sphere,ica_compvars,~,~,~,ica_act_r] = runica(deegr','pca',ncomps2keep);
%     ica_winv = pinv(ica_weights*ica_sphere); % pseudo inverse to get forward (?) weights for topo-plotting
%     ica_ncomps = length(ica_compvars);
%     ica_act = reshape(ica_act_r, ica_ncomps, nsamps, ntrials);
% 
%     % run GED 
%     %!
%     
%     % ?? across trials would be worth trying too, maybe in a different script though
% 
%     %% LINEAR MODEL - UNIVARIATE W/ EPOCH TIMEPOINTS AS OBSERVATIONS
% 
%     % select component type to use:
% %     comptype = 'PCA';
%     comptype = 'ICA';
% %     comptype = 'GED';
%     %
% 
%     % get components to use for the model
%     switch comptype
%         case 'PCA'
%             coef2use   = pca_coef;  
%             scores2use = pca_scores;
%         case 'ICA'
%             coef2use   = ica_winv;
%             scores2use = ica_act;
%         case 'GED'
%             %!!!!!
%     end
%     ncomps = size(coef2use,2);
% 
%     % loop through components and plot topographies 
%     topo = coef2use;
%     figure('name','component topographies');
%     for k = 1:ncomps
%         subplot(5,5,k);
%         rs_topoplot(topo(:,k), chanlocs, 'colormap',colormap2); 
%         title(num2str(k))
% %         c = 5; clim([-c,c])
%     end
% 
%     % handle data
%     dlfp2      = squeeze(dlfp);
%     dlfp2_r    = reshape(dlfp2,nsamps*ntrials,1);
%     scores2use = permute(scores2use,[3,2,1]);
% 
%     % get electrode locations to use as a grouping variable
%     uni = unique(tbl_depths2.dist_angle);
% 
%     % repeat table once per trial/timepoint pair
%     inds = repmat((1:ntrials)', nsamps,1);
%     tbl_depths3 = tbl_depths2(inds,:);
% 
%     % define models to run
%     ctemp = [];
%     ctemp.Y = 'eeg';
%     ctemp.X_tar = 'lfp'; % target predictor
%     ctemp.X_ints = {'Depth','DepthRelTOA'}; % interaction terms for X_tar
% %     ctemp.X_ints = {'Depth','DepthRelTOA', 'DepthRelDura'};
%     ctemp.G = {'dist_angle'}; % grouping variables
% %     ctemp.G = {'dist_angle', 'Depth_2split'}; % ! Depth-based median split or N split
%     models  = genLMEformulas(ctemp);
%     nmodels = length(models)
% % ? models to add:
% %   ? another way to divide electrodes, e.g. forward back location, medial lateral location
% %   ? median- or n-split by Depth, DepthRelTOA and DepthRelDura
% 
%     % components as separate models w/ timepoints as observations
%     mdls = cell(ncomps,1);
%     for c = 1:ncomps
% %     c = 3;
% 
%         scores2use_temp = squeeze(scores2use(:,:,c))';
%         scores2use_temp = reshape(scores2use_temp, nsamps*ntrials,1);
% 
%             tbl = tbl_depths3;
%             tbl.eeg = scores2use_temp;
%             tbl.lfp = dlfp2_r;
% 
%             % loop through models
%             mdls_temp = cell(nmodels,1);
%             bic_temp = nan(nmodels,1);
%             parfor md = 1:nmodels
%                 mdls_temp{md} = fitlme(tbl, models{md} ); % intercept only
%                 bic_temp(md) = mdls_temp{md}.ModelCriterion.BIC;
%             end
%             [~, best_model] = min(bic_temp); 
%             mdls(c,1)  = mdls_temp(best_model);
%             
%             % !!! NEED TO STORE MORE THAN JUST BEST MODEL, SINCE several models have fairly similar BIC
%             % !!! IF FIT IS SINGULAR, BIC CANNOT BE USED AS POSTERIOR CAN'T BE ESTIMATED BY NORMAL DISTRIBUTIONS
%             
%           c
%     end % COMP LOOP
% 
%     %% PLOT
%     close all
% 
%     % plot component topographies with best fitting models titles
%     topo = coef2use;
%     figure('name','component topographies');
%     for c = 1:ncomps
%         subplot(5,5,c);
%         rs_topoplot(topo(:,c), chanlocs, 'colormap',colormap2); 
%         m = mdls{c};
%         
%         title(num2str( m.Formula.char ))
% %         title(num2str( m.ModelCriterion.BIC ))
% %         title(num2str( m.Rsquared.adjusted ))
% %         title(num2str( max(abs(m.Coefficients.Estimate(contains(m.CoefficientNames,'lfp')))) ))
% %         % ? are amplitudes of each score variable the same? affects interpretation!!!
%         % THEY ARE NOT THE SAME, NORMALISE THEM FOR COMPARISON? 
% 
% %         cl = 5; clim([-cl,cl])
%     end
% 
%     % extract info
%     rsqa = nan(ncomps,1);
%     bic  = nan(ncomps,1);
%     bestmodels = cell(ncomps,1);
%     fixed_est = cell(ncomps,1);
%     fixed_ci  = cell(ncomps,1);
%     fixed_names = cell(ncomps,1);
%     rand_est   = cell(ncomps,1);
% %     rand_names = cell(ncomps,1); 
%     for c = 1:ncomps
%         m = mdls{c};
%         bestmodels{c} = m.Formula.char;
%         rsqa(c) = m.Rsquared.adjusted;
%         bic(c)  = m.ModelCriterion.BIC;
%         fixed_est{c} = m.Coefficients.Estimate;
%         fixed_names{c} = m.CoefficientNames;
%         fixed_ci{c}    = m.coefCI;
%         try 
%             rand_est{c} = m.randomEffects;
% %             rand_names{c} = m. !!
%         end
% 
%     end
% 
%     % bar plot of R-squared
%     figure; bar(rsqa); title 'R-squared (adjusted)'
%     
%     % best-fitting model component coefficients
%     figure('name','best-fitting model component coefficients');
%     for c = 1:ncomps
%         subplot(5,5,c);
%         bar(fixed_est{c}); hold on;
%         errorbar(fixed_est{c},[fixed_ci{c}(:,1), fixed_ci{c}(:,2)]) % !CHECK THIS WORKS PROPERLY FOR MULTIPLE PREDICTORS
%         ax = gca; ax.XTickLabels = fixed_names{c};                  %!!!? why two sets of error bars?
% %         cl = 0.5 + ceil(max(abs([fixed_est{:}]))); ylim([-cl,cl]); 
%     end
%     % ? WHY are models often fitting with a very precisely zero fixed or random effect? should have models with no fixed effects or even intercept? but why random does that too sometimes...
% 
%     % best-fitting model random effects
%     figure('name','best-fitting model component coefficients');
%     for c = 1:ncomps
%         subplot(5,5,c);
%         bar(rand_est{c});
% %         cl = 0.5 + ceil(max(abs([fixed_est{:}]))); ylim([-cl,cl])
%         if isempty(rand_est{c})
%             title '* NO RAND EFFECTS *'
%         else
% %             title 
%         end
%     end
% 
%     % compute metric of "widespread-ness" for each component topography 
%     wideness = nan(ncomps,1);
%     wideness_fisher = nan(ncomps,1);
%     for c = 1:ncomps
%         topo = coef2use(:,c);
%         cov = topo * topo'; % compute covariance matrix
%         % normalise to population pearson coefficient
%         for row = 1:nchans 
%             for col = 1:nchans
%                 cov(row,col) = cov(row,col) ./ (sqrt(cov(row,row))*sqrt(cov(col,col)));
%             end
%         end
%         % sum of absolute correlations, excluding self-correlation, divided by perfect correlation of all electrodes
%         wideness(c) = (sum(abs(cov),'all')-nchans) / nchans^2; 
%         wideness_fisher(c) = atanh(wideness(c)); % Fisher transformation for linear/normal fitting
%     end
%     
%     % sort by widespreadness metric and plot component topographies
%     [wideness_sorted, sortinds_w] = sort(wideness,'descend');
%     wideness_fisher_sorted = wideness_fisher(sortinds_w);
%     coef2use_sorted_w = coef2use(:,sortinds_w);
%     figure('name','component topographies sorted by widespreadness metric');
%     for c = 1:ncomps
%         subplot(5,5,c);
%         rs_topoplot(coef2use_sorted_w(:,c), chanlocs, 'colormap',colormap2); 
%         title(num2str( wideness_sorted(c) ))
% %         cl = 5; clim([-cl,cl])
%     end
%     
%     % plot best-fit model R^2 against component variance and widespreadness
%     figure('name','R-squared (adjusted) vs comp variance & widespreadness');
%     % compute component variance
%     scores2use_r = reshape(scores2use,ntrials*nsamps,ncomps);
% %     compvarmet = sum(abs(scores2use_r))';
% %     compvarmet = mean(abs(scores2use_r))';
% %     compvarmet = std(abs(scores2use_r))';
%     compvarmet = nan(ncomps,1);
%     for c = 1:ncomps
%         [~, compvarmet(c)] = compvar( deegr', scores2use_r', coef2use, c);
%     end
%     compvarmet = compvarmet * -1;
%     100 * compvarmet / sum(compvarmet)
%     %
%     subplot(1,2,1); scatter(compvarmet,rsqa); axis square; lsline; 
%                     xlabel 'component "variance" metric'; ylabel 'Rsquared (adjusted)'; 
%                     title(sprintf('R = %.2f',corr(compvarmet,rsqa)))
%     subplot(1,2,2); scatter(wideness_fisher,rsqa); axis square; lsline; 
%                     xlabel 'wideness score (Fisher transformed)'; ylabel 'Rsquared (adjusted)'; 
%                     title(sprintf('R = %.2f',corr(wideness_fisher,rsqa)))
%     
%     % ? of course some small variance, localised components could be well predicted by LFP if nearby to recording site
%     %       can test by plotting Rsqa against some metric of topo variance scaled by proximity to recording site electrode (check notes)
% 
% 
%  
% %% END LOOP THROUGH CONDS
% end


%% END
% cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


