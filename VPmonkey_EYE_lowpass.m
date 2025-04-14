% 
% 
% 
%   apply lowpass filter to EYE data 
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

s.dsf = 2;

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
files = files(startsWith(files,'merged_EYE'));
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

    %% switch to eeglab
    evalc 'addpath(genpath(pathlab));';
    evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''fieldtrip-20210311'' ]))'; % remove fieldtrip to tidy things up a bit
    evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''Biosig3.7.9'' ]))'; % remove biosig to tidy things up a bit
    evalc 'rmpath(genpath(pathlw));';

    %% load file 
    EYE = pop_loadset(files{f}); 

    %% LOWPASS FILTER
    chan_pup = find(strcmp({EYE.chanlocs.labels},'eye_PupilDiameter'));
    PUP = pop_select(EYE,'channel', chan_pup);
    PUP = pop_eegfiltnew( PUP, 'hicutoff', s.lowpass);
    EYE.data(chan_pup,:,:) = PUP.data(1,:,:);

    %% SAVE FILTERED DATA
    saveName = ['LP_' num2str(s.lowpass) ' ' files{f}];
    EYE = pop_saveset(EYE, 'filepath', homedir  ...
        ,'filename', saveName ); 

    %% export to lw
    cd([homedir filesep 'lw'])

    % downsample for lighter lw storage
    EYE = pop_rs_downsample(EYE,s.dsf);
    
    % crop for lw
    EYE = pop_select(EYE, 'time', s.xlims_export );

    % export
    cd([s.savePath filesep 'lw'])
    saveName = strrep(saveName,'.set','');
    rs_convert_lab2lw_V1( EYE, saveName, [] );
    fprintf('... saved file %s: %s\n', saveName, [cd filesep saveName])  

    %% split by condition & export to lw
    ents = EYE.event;
    conds = unique({ents.type});
    for c = 1:length(conds)

        % extract this condition
        cond = conds{c}

        EEG_cond = pop_select(EYE,'trial', find(strcmp({ents.type},cond)) );
        
        % replace all NaNs with zeros for easier lw plotting
        EEG_cond.data(isnan(EEG_cond.data)) = 0;

        % export to lw
        rs_convert_lab2lw_V1( EEG_cond ...
            , [ 'ep ' cond ' ' saveName ], [] );

    end

    %% switch to letswave
    evalc 'rmpath(genpath(pathlab));';
    evalc 'addpath(genpath(pathlw));';

    %% compute averages of each condition after ICA-cleaning
    for c = 1:length(conds)

        cond = conds{c}

        % add channel locations
        filename = [ s.savePath filesep 'lw' filesep 'ep ' cond ' ' EYE.filename(1:end-4) '.lw6'];
        option=struct('filename', filename);
        lwdata = FLW_load.get_lwdata(option);
        option=struct('filepath','/home/rick/Dropbox/Somervail_Richard/toolbox/letswave7-master/res/electrodes/spherical_locations/monkeyEEG_LW.locs','suffix','','is_save',1);
        lwdata= FLW_electrode_location_assign.get_lwdata(lwdata,option);

        % take average
        option=struct('filename', filename);
        lwdata = FLW_load.get_lwdata(option);
        option=struct('operation','average','suffix','avg','is_save',1);
        lwdata = FLW_average_epochs.get_lwdata(lwdata,option);

    end
    

end % end session loop

cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

