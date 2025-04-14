% 
% 
% 
%   apply lowpass filter to LFP data before plotting & CCA
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

s.lowpass = 30; % to match EEG low pass filter

% s.dsf = 2;

% % channel locations for these datasets
% s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
% s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
% s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
% chanlocs = pop_readlocs(s.chanlocs); 
% chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
% chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
% load(s.chanlocs_lw);

%% select files 
cd(homedir)
files = dir; files = {files.name}'; 

%%%%%%% SELECT FILE TYPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
files = files(startsWith(files,'merged_LFP'));
files = files(endsWith(files,'.set'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

files(contains(files,'all')) = []; % exclude files which contain data merged across sessions
% files = files(endsWith(files,'M.set'));
% files = files(endsWith(files,'T.set'));
if isempty(files), error 'RS error: no files found'; end

%
sel = listdlg('ListString', files, 'SelectionMode','multiple', 'ListSize',[400,600]);
% sel = 2; warning( '! BYPASSING SESSION SELECTION' )
files = files(sel);

nfiles = length(files);

%% loop through sessions
for f = 1 :nfiles

    cd(homedir)

    %% load file 
    LFP = pop_loadset(files{f}); 

    %% LOWPASS FILTER
    LFP = pop_eegfiltnew( LFP, 'hicutoff', s.lowpass);

    %% SAVE FILTERED DATA
    saveName = ['LP_' num2str(s.lowpass) ' ' files{f}];
    LFP = pop_saveset(LFP, 'filepath', homedir  ...
        ,'filename', saveName ); 

end % end session loop

cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

