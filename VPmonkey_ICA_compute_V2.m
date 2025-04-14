%
%  compute ICA (but no IC rejection yet)
%
%
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
addpath([getRoot '/VPmonkey/scripts'])

% try 
%     eeglab
% catch
%     addlab
% end

%% SETTINGS
s = [];
s.savePath =  homedir; % mkdir(s.savePath)

% ica
s.useICLabel   =  false;
s.useBlinkCorr =  false;
s.ica_extended =  true;
% s.ica_extended = false; % extended=1 is generally recommended for better isolating subgaussian ICs like line-noise and slow drifts
% s.ica_stop = 1E-6; % default
s.ica_stop = 1E-7; % longer for better decomposition

% % channel locations for these datasets
% s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
% s.chanlay  = '/home/rick/neuro/iannettilab/Projects/VP monkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
% s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
% chanlocs = pop_readlocs(s.chanlocs); 
% chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
% chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
% load(s.chanlocs_lw);


%% select files to ICA
files = dir; files = {files.name}'; 
files = files(startsWith(files,'noICA') & endsWith(files,'.set') & contains(files,'merged_EEG') );
files(contains(files,'all')) = []; % exclude files which contain data merged across sessions
if isempty(files), error 'RS error: no files found'; end
sel = listdlg('ListString', files, 'SelectionMode','multiple', 'ListSize',[400,600]);
% sel = 2; warning( '! BYPASSING SESSION SELECTION' )
files = files(sel);
nfiles = length(files);

%% loop through sessions
for f = 1 :nfiles

    %% load file 
    EEG = pop_loadset(files{f}); 
    
    %% determine which channels to use for ICA based on the reference used, then only use those channels for ICA
    numChansInterpolated = EEG.etc.numChansInterpolated;
    chans = {EEG.chanlocs.labels};
    types = {EEG.chanlocs.type};
     
    if ~any(startsWith( strsplit(EEG.ref)  ,'EXG')) % looks like CAR or ring reference, so keep EXG chans if they are present 
%         EEG_otherchans = pop_select(EEG, 'channel', find(~strcmp({EEG.chanlocs.type},'EEG') & ~startsWith({EEG.chanlocs.labels},'EXG')) );
%         EEG_4ICA       = pop_select(EEG, 'channel', find( strcmp({EEG.chanlocs.type},'EEG') |  startsWith({EEG.chanlocs.labels},'EXG')) );
        chanindy_4ica = find( strcmp(types,'EEG') | ( startsWith(chans,'EXG') & ~strcmp(chans,'EXG3') )  ); % never use headpost in ICA though
        numIC = length(chanindy_4ica) - numChansInterpolated - 1;  % remove one additional IC for the loss of rank from the CAR/ring ref
    elseif any(startsWith( strsplit(EEG.ref)  ,'EXG')) % if EXG channels were used as reference, then remove them from ICA
%         EEG_otherchans = pop_select(EEG, 'channel', find(~strcmp({EEG.chanlocs.type},'EEG'))  );
%         EEG_4ICA       = pop_select(EEG, 'channel', find( strcmp({EEG.chanlocs.type},'EEG'))  );
        chanindy_4ica = find( strcmp(types,'EEG') & ~startsWith(chans,'EXG')  ); % don't include ears or headpost
        numIC = length(chanindy_4ica) - numChansInterpolated ;  % number of ICs should just be number of EEG channels minus number of interpolated channels
    else 
        error 'unclear what reference was used ... investigate and fix conditions'
    end

    %% compute ICA
    EEG  = pop_runica(EEG, 'icatype','runica', 'stop', s.ica_stop, ...
        'pca', numIC, 'extended', s.ica_extended, 'chanind', chanindy_4ica );

    %% use IClabel % ? of course this may not be so accurate when used on a different species to what it was trained on 
    if s.useICLabel
        EEG = iclabel(EEG);
    end

    %% save 
    EEG.etc.preproc_ica = s; % save ICA settings here

    % save .set
    if      startsWith(EEG.setname, 'noICA')
        saveName = [ strrep( EEG.setname, 'noICA', 'cpICA') ];

        % delete the noICA version of file (this is already exported to lw anyway)
        delete(files{f}); % delete set
        delete(strrep(files{f},'.set','.fdt')); % delete fdt
    elseif  startsWith(EEG.setname, 'cpICA')
        saveName = EEG.setname;
    end
    EEG.setname  =  saveName;
    EEG.filename = [saveName '.set'];
    pop_saveset(EEG, 'filepath', s.savePath, 'filename', EEG.filename);
    fprintf('completed session "%s" (file %d/%d)\n\n\n', files{f}, f, nfiles)

end % end session loop

cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

