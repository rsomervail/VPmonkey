% 
% 
%   
%   apply GED temporal filter to remove signal dropouts then lowpass filter to EYE data 
%   before subsequent analysis
% 
%       UNFINISHED - MIGHT SCRAP AND JUST DO ONLINE AFTER ALL
%
%       - Richard Somervail, 2024
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
addpath(genpath([getRoot '/VPmonkey/scripts']))

pathlw  = [getRoot 'toolbox' filesep 'letswave7-master'];
pathlab = [getRoot 'toolbox' filesep 'eeglab2022.1'];

%% switch to eeglab
evalc 'addpath(genpath(pathlab));';
evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''fieldtrip-20210311'' ]))'; % remove fieldtrip to tidy things up a bit
evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''Biosig3.7.9'' ]))'; % remove biosig to tidy things up a bit
evalc 'rmpath(genpath(pathlw));';

%% SETTINGS
s = [];
s.savePath =  homedir; % mkdir(s.savePath)

s.xlims_export  = [-0.5 2];

s.plotYlims2 = [-10 10]; % plot limits for ft_databrowser (rejecting epochs)

s.lowpass = 5; % could try even lower like 2 Hz but might be a bit too extreme

s.dsf = 2; % only for lw export

% % channel locations for these datasets
% s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
% s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
% s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
% chanlocs = pop_readlocs(s.chanlocs); 
% chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
% chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
% load(s.chanlocs_lw);



subs = {'SubM','SubT'};

for sub = subs

    %% load EYE data
    cfg = struct; cfg.exportlw = false;
    cfg.sub = sub{1};
    cfg.filetypes = {'EYE'};
    cfg.LFP_include = {'LP_30'};
    cfg.filt_eye = false; % preproc script needs all sessions
    cfg.EYE_exclude = {'LP_','GED_TEMPO'}; % exclude any previously low-pass filtered or GED temporal cleaned data
    cfg.byElectrode = false;
    cfg.average   = false;
    cfg.mergesesh = false;
    cfg.zscore = {'EYE'};  % z-scoring is necessary here for appropriate merging, so this detail is listed in output filename
%     cfg.zscore_win.EYE = [-1 3]; % ? use whole epoch for now
    cfg.zscore_cond  = false; % z-scoring done across all conditions at this stage, it can always be redone per condition if needed
    data = VPmonkey_mergeSesh(cfg);
    
    nfiles = length(data);
    
    %% CONCAT PUPIL DATA ACROSS ALL SESSIONS
    d = double(cat(3,data.EYE.data));
    
    
    
    % !! UNFINISHED FROM HERE
    %  ? probs can take inverse filter and apply to each session individually?
    %!! when exporting, I need to rename each file appropriately and specify that it is PUP only now

end

% OLD CODE BELOW

%% loop through sessions
for f = 1 :nfiles

    cd(homedir)

    %% load file 
    EYE = data.EYE(f)

    %% GED TEMPORAL FILTER
    chan_pup = find(strcmp({EYE.chanlocs.labels},'eye_PupilDiameter'));
    PUP = pop_select(EYE,'channel', chan_pup);
    PUP = pop_rs_removeNonZeroEvents(PUP);

    % loop through conditions
    conds_all = {PUP.event.type}; conds = unique(conds_all); nconds = length(conds);
    dpup = squeeze(double(PUP.data));
    for cond = 1:nconds

        % get trials belonging to this condition
        inds = find(strcmp(conds_all,conds{cond}));
        dpup_cond = dpup(:,inds);
    
        % compute temporal GED filter
        mu = mean(dpup_cond,2);
        dpup_cond_c = dpup_cond - mu;
        [evals, inverse, forward, scores] = rs_ged(data, cfg)

        true

    end

    % !! PLOT BEFORE/AFTER FOR THIS SESSION AND SAVE FIGURE

    % STORE CLEANED VERSION AS OUTPUT
    EYE.data(chan_pup,:,:) = PUP.data(1,:,:);


    


    %% LOWPASS FILTER
    chan_pup = find(strcmp({EYE.chanlocs.labels},'eye_PupilDiameter'));
    PUP = pop_select(EYE,'channel', chan_pup);
    PUP = pop_eegfiltnew( PUP, 'hicutoff', s.lowpass);
    EYE.data(chan_pup,:,:) = PUP.data(1,:,:);

    %% SAVE FILTERED DATA
    saveName = ['LP_' num2str(s.lowpass) 'GED_TEMPO ' files{f}];
    EYE = pop_saveset(EYE, 'filepath', homedir  ...
        ,'filename', saveName ); 

end % end session loop

cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

