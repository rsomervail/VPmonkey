% 
% 
%   Batch script to use VPmonkey_mergeSesh function to load/export data in several formats
% 
% 
%   
%%
clc
clearvars
close all

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

homedir = ['/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_' subfold];
cd(homedir);

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

subs = {'SubM','SubT'};

%% RUN BATCH
tic
for sb = 1:length(subs)
    sub = subs{sb};

%% NO TRIAL REJECTION
% 
% % GOOD EEG
% cfg = struct; cfg.exportlw = true;
% cfg.sub = sub;
% cfg.filetypes = {'EEG', 'EYE'};
% cfg.byElectrode = false;
% cfg.average   = true;
% cfg.mergesesh = true;
% cfg.filt_preproc = 'allGood'; % 'allSesh'
% cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
% cfg.EEG_exclude = {'noICA'};
% data_sesh = VPmonkey_mergeSesh(cfg);
% % % ICA ONLY
% % cfg = struct; cfg.exportlw = true; 
% % cfg.sub = sub;
% % cfg.filetypes = {'EEG', 'EYE'};
% % cfg.byElectrode = false;
% % cfg.average   = true;
% % cfg.mergesesh = true;
% % cfg.filt_preproc = 'allGood'; 
% % cfg.EEG_include = {'yesICA'};
% % cfg.EEG_exclude = {'noICA','LP_30','ged_horiz'};
% % cfg.EYE_include = {'LP_5'};
% % data_sesh = VPmonkey_mergeSesh(cfg);
% 
% % ALL EEG 
% cfg = struct; cfg.exportlw = true;
% cfg.sub = sub;
% cfg.filetypes = {'EEG', 'EYE'};
% cfg.byElectrode = false;
% cfg.average   = true;
% cfg.mergesesh = true;
% cfg.filt_preproc = 'allSesh';
% cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
% cfg.EEG_exclude = {'noICA'};
% data_sesh = VPmonkey_mergeSesh(cfg);
% % % ICA ONLY
% % cfg = struct; cfg.exportlw = true;
% % cfg.sub = sub;
% % cfg.filetypes = {'EEG', 'EYE'};
% % cfg.byElectrode = false;
% % cfg.average   = true;
% % cfg.mergesesh = true;
% % cfg.filt_preproc = 'allSesh';
% % cfg.EEG_include = {'yesICA'};
% % cfg.EEG_exclude = {'LP_30','ged_horiz'};
% % cfg.EYE_include = {'LP_5'};
% % data_sesh = VPmonkey_mergeSesh(cfg);
% 
% % by electrode (LFP, MUA etc)
% cfg = struct; cfg.exportlw = true;
% cfg.sub = sub;
% % cfg.filetypes = {'MUA','LFP'};
% cfg.filetypes = {'LFP'};
% cfg.byElectrode = true;
% cfg.average   = true;
% cfg.mergesesh = true;
% cfg.zscore = cfg.filetypes; 
% cfg.zscore_win = [-0.2 0.6];
% cfg.zscore_cond = true;
% data_elec = VPmonkey_mergeSesh(cfg);

%% WITH TRIAL REJECTION

ar = struct;
ar.method = 'time';
ar.metric = 'median';
ar.thresh = 3;
ar.timeprop = 0.1; % ? too strict?
ar.chanprop = 0.1; % ? too strict?

% % GOOD EEG
% cfg = struct; cfg.exportlw = true;
% cfg.autoar = ar;
% cfg.sub = sub;
% cfg.filetypes = {'EEG', 'EYE'};
% cfg.byElectrode = false;
% cfg.average   = true;
% cfg.mergesesh = true;
% cfg.filt_preproc = 'allGood'; 
% cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
% cfg.EEG_exclude = {'noICA'};
% data_sesh = VPmonkey_mergeSesh(cfg);
% % % ICA ONLY
% % cfg = struct; cfg.exportlw = true;
% % cfg.sub = sub;
% % cfg.filetypes = {'EEG', 'EYE'};
% % cfg.byElectrode = false;
% % cfg.average   = true;
% % cfg.mergesesh = true;
% % cfg.filt_preproc = 'allGood'; 
% % cfg.EEG_include = {'yesICA'};
% % cfg.EEG_exclude = {'noICA','LP_30','ged_horiz'};
% % data_sesh = VPmonkey_mergeSesh(cfg);

% % ALL EEG 
% cfg = struct; cfg.exportlw = true;
% cfg.autoar = ar;
% cfg.sub = sub;
% cfg.filetypes = {'EEG', 'EYE'};
% cfg.byElectrode = false;
% cfg.average   = true;
% cfg.mergesesh = true;
% cfg.filt_preproc = 'allSesh'; 
% cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
% cfg.EEG_exclude = {'noICA'};
% cfg.EYE_include = {'LP_5'};
% data_sesh = VPmonkey_mergeSesh(cfg);
% % ICA ONLY
% cfg = struct; cfg.exportlw = true;
% cfg.sub = sub;
% cfg.filetypes = {'EEG', 'EYE'};
% cfg.byElectrode = false;
% cfg.average   = true;
% cfg.mergesesh = true;
% cfg.filt_preproc = 'allSesh'; 
% cfg.EEG_include = {'yesICA'};
% cfg.EEG_exclude = {'noICA','LP_30','ged_horiz'};
% data_sesh = VPmonkey_mergeSesh(cfg);

% by electrode (LFP, MUA etc)
cfg = struct; cfg.exportlw = true;
cfg.autoar = ar;
cfg.sub = sub;
% cfg.filetypes = {'MUA','LFP'};
cfg.filetypes = {'LFP'};
cfg.byElectrode = true;
cfg.average   = true;
cfg.mergesesh = true;
cfg.zscore = cfg.filetypes; 
cfg.zscore_win = [-0.2 0.6];
cfg.zscore_cond = true;
cfg.LFP_include = {'LP_30'};
data_elec = VPmonkey_mergeSesh(cfg);

end
toc
