% 
% 
%       ! important to not throw away uncertainty in level 1 parameter estimates here..
%         maybe can do this by sampling from parameter estimates given level 1 SE?
% 
% 
%%
clearvars

sub = 'SubM';

load([getRoot '/VPmonkey/paper/results/'  sub '_results_LM_LFP_EEG_LEVEL1.mat'])

% !! REMOVE AFTER RE-RUNNING LEVEL 1 SCRIPT ONCE
r = results; clear results 
r.conds = {'VIS', 'SOM', 'AUD'}; nconds = length(r.conds);
r.LFP_est = r.est;
r.LFP_SE  = r.SE;
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% 
cond = 3;

rsq2plot = sqz(r.rsq(cond,:,:,:));

figure; 
subplot(2,1,1); plot( squeeze(mean(rsq2plot,'omitnan'))' ); title(r.conds{cond})
subplot(2,1,2); plot( squeeze(mean(rsq2plot,1:2,'omitnan'))' );

%% FIT MODELS TO LEVEL 1 LFP COEFFICIENT ESTIMATES

t = 157; % ! Rsquared peak, no idea what corresponds to because I don't have the time vector yet

% format variables for model
LFP_est = squeeze( r.LFP_est(cond,:,:,t));
LFP_SE  = squeeze(  r.LFP_SE(cond,:,:,t));
[nelecs,ncomps] = size(LFP_est);
LFP_est = reshape(LFP_est,nelecs*ncomps,1);
LFP_SE  = reshape(LFP_SE,nelecs*ncomps,1);
weights = 1./LFP_SE;
%
elec_depths  = r.tbl_depths.DepthRelTOA;
elec_depths  = repmat(elec_depths, ncomps,1);
%
comp_ID      = repmat((1:ncomps)', 1, nelecs)';
comp_ID      = reshape(comp_ID,nelecs*ncomps,1);
%
comp_wideness = atanh(r.wideness);
comp_wideness = repmat(comp_wideness, 1, nelecs)';
comp_wideness = reshape(comp_wideness,nelecs*ncomps,1);

tbl = table(LFP_est,comp_wideness, comp_ID, elec_depths);
tbl.comp_ID = categorical(tbl.comp_ID);

% ! CURRENTLY REMOVING NANS FROM TIMEPOINT / COMP PAIRS WHERE THE BEST FIT DID NOT INCLUDE LFP PREDICTOR ?? I THINK
% ? ANYWAY FOR NOW REMOVING THESE, LATER ON COULD JUST INCLUDE THEM SINCE THEY WOULD BE LOW-WEIGHTED 
badrows = isnan(tbl.LFP_est);
tbl(badrows,:)   = [];
weights(badrows) = [];

% define model formulas
% form = 'LFP_est ~ comp_wideness * elec_depths';
% form = 'LFP_est ~ -1 + comp_ID';
form = 'LFP_est ~ -1 + comp_ID + (elec_depths - 1 | comp_ID)';

% ? maybe this is too complicated, and need to do electrode depths and component analysis as separate levels?
% e.g. one model per component 

% fit model
m = fitlme(tbl,form,'Weights',weights, 'DummyVarCoding', 'full')

% ? maybe should have random slope per component to deal with discereteness problem
%  & also probable issue with arbitrary signs!!!
% ? deals with issue of discreteness in component metrics
%     - can then do 3rd level analysis where I plot component metrics against mean slopes??
% ? similar problem also with recording location...

% ? should still eventually have way to explicit observation uncertainty, i.e. an explicit multi-level model


%% PLOT 

% plot all coefficients 
figure; errorbar( m.Coefficients.Estimate, m.Coefficients.SE, 'o'); refline(0,0); 
xticks(1:length(m.CoefficientNames)); xticklabels(strrep(m.CoefficientNames,'_','-'))

% plot mean model coefficient per component against component metrics
figure; scatter( r.wideness, m.Coefficients.Estimate); xlabel 'comp wideness'; ylabel 'mean LFP-est'; lsline
% figure; errorbar( r.wideness, m.Coefficients.Estimate, m.Coefficients.SE,'o'); xlabel 'comp wideness'; ylabel 'mean(abs(LFP-est)'; 


% figure; scatter( tbl.comp_wideness, tbl.LFP_est ); refline(0,0);
% % ? would be nice to include weights here somehow, opacity? error bars?





