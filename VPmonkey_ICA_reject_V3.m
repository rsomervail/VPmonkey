%
%  reject ICA components (using fancy new GUI plugin I made)
%
%       - Richard Somervail, 2023
%%
clc
clearvars
close all

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

homedir = ['/media/rick/Rick_LabData_3/neuro/iannettilab_ongoing/VPmonkey/data/clean_' subfold];
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

s.freqs = [1,100];

s.dsf = 2;

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

%% select files to ICA
files = dir; files = {files.name}'; 
files = files(startsWith(files,'cpICA') & endsWith(files,'.set') & contains(files,'merged_EEG') );
files(contains(files,'all')) = []; % exclude files which contain data merged across sessions
files = [ files(endsWith(files,'M.set')) , files(endsWith(files,'T.set'))  ];
if isempty(files), error 'RS error: no files found'; end
sel = listdlg('ListString', files, 'SelectionMode','multiple', 'ListSize',[400,600]);
% sel = 2; warning( '! BYPASSING SESSION SELECTION' )
files = files(sel);
nfiles = length(files);

%% loop through sessions
for f = 1 :nfiles

    %% switch to eeglab
    evalc 'addpath(genpath(pathlab));';
    evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''fieldtrip-20210311'' ]))'; % remove fieldtrip to tidy things up a bit
    evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''Biosig3.7.9'' ]))'; % remove biosig to tidy things up a bit
    evalc 'rmpath(genpath(pathlw));';

    %% load file 
    cd(homedir)
    EEG = pop_loadset(files{f}); 
    numIC = size(EEG.icaweights,1); % get number of ICA componentts in this dataset

    if isfield( EEG.etc, 'ic_classification') %  IClabel is not helpful here.
        EEG.etc = rmfield( EEG.etc, 'ic_classification');
    end

    %% check for previously-cleaned version of file 
    file_cleaned = strrep(files{f},'cpICA','yesICA');
    if exist( file_cleaned,'file')
        info = pop_loadset( 'filename', file_cleaned,'loadmode','info' );
        comps2reject = info.etc.preproc.ica.rejIC;
    else
        comps2reject = 1;
    end
    
    %% reject ICs ----------------------------------------------------------------------------------------------------------
    close all;
    %
    ctemp = [];
    ctemp.letsplot.channel = 'CZ';
    ctemp.letsplot.xlim    = [-0.2, 0.6];
%     ctemp.letsplot.ylim    = [-20, 20];
    ctemp.freqs         = s.freqs;
    ctemp.comps2reject  = comps2reject; 
    ctemp.arrangefigs   = false;
    %
    EEG = pop_rs_ica_reject(EEG,ctemp);


%     % BYPASS GUI AND JUST APPLY CLEANING WITH CURRENTLY SELECTED COMPONENTS 
%     if exist( file_cleaned,'file')  % (FOR EXPORTING TO LW)
%         temp = update_rs_ica_applyCleaning(EEG, comps2reject );
%         EEG = temp(2);
%         EEG.etc.preproc.ica.rejIC = comps2reject;
%         EEG.etc.preproc.ica.numIC = numIC;
%     end

    % ----------------------------------------------------------------------------------------------------------------------
 
    %% IF ANY REJECTION WAS ACTUALLY APPLIED
    if ~isnan(EEG.etc.preproc.ica.rejIC)

        %% save .set 
        saveName = [ strrep( EEG.filename(1:end-4), 'cpICA', 'yesICA') ];
        EEG.setname  = saveName;
        EEG.filename = [saveName '.set'];
        pop_saveset(EEG, 'filepath', s.savePath, 'filename', EEG.filename);

        %% export to lw
    
        % downsample for lighter lw storage
        EEG = pop_rs_downsample(EEG,s.dsf);
        
        % crop for lw
        EEG = pop_select(EEG, 'time', s.xlims_export );
    
        % export
        cd([s.savePath filesep 'lw'])
        rs_convert_lab2lw_V1( EEG, saveName, chanlocs_lw );
        fprintf('... saved file %s: %s\n', saveName, [cd filesep saveName])  

        %% split by condition & export to lw
        ents = EEG.event;
        conds = unique({ents.type});
        for c = 1:length(conds)
    
            % extract this condition
            cond = conds{c}

            EEG_cond = pop_select(EEG,'trial', find(strcmp({ents.type},cond)) );
            
            % replace all NaNs with zeros for easier lw plotting
            EEG_cond.data(isnan(EEG_cond.data)) = 0;
    
            % export to lw
            saveName = [ 'ep ' cond ' ' EEG.filename(1:end-4) ];
            rs_convert_lab2lw_V1( EEG_cond ...
                , saveName, [] );
    
        end

        %% switch to letswave
        evalc 'rmpath(genpath(pathlab));';
        evalc 'addpath(genpath(pathlw));';

        %% compute averages of each condition after ICA-cleaning
        for c = 1:length(conds)

            cond = conds{c}

            % add channel locations
            filename = [ s.savePath filesep 'lw' filesep 'ep ' cond ' ' EEG.filename(1:end-4) '.lw6'];
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



    %% if no rejection was applied
    else
        fprintf('no ICs were rejected ... cancelling\n')
    end
    

end % end session loop

cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

