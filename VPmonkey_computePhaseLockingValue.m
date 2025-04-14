%
%    run phase-locking analysis on each condition merged across all files
% 
% 
%%
clc
clearvars
close all

subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
% subfold = 'main_example';

homedir = ['/media/rick/Rick_LabData/neuro/iannettilab_ongoing/VPmonkey/data/clean_' subfold]; % external HDD
% homedir = ['/home/rick/neuro/iannettilab/Projects/VP monkey/data/clean_' subfold]; % internal SSD
cd(homedir);

addpath('/home/rick/neuro/iannettilab/Projects/VP monkey/Giac_ToolBox')
addpath('/home/rick/neuro/iannettilab/Projects/VP monkey/scripts')

% try 
%     eeglab
% catch
%     addlab
% end

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/neuro/iannettilab/Projects/VP monkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

%% SETTINGS
s = [];
% s.savePath =  ['/media/rick/Rick_LabData/neuro/iannettilab_ongoing/VPmonkey/data/clean_' subfold]; % mkdir(s.savePath)
s.savePath =  [homedir]; % mkdir(s.savePath)

s.conds = {'VIS', 'SOM', 'AUD'};
s.epwin = [ -2  2 ]; % window of interest

s.xlims_export = [-1 2];

% % channel locations for these datasets
% s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
% s.chanlay  = '/home/rick/neuro/iannettilab/Projects/VP monkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
% s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
% chanlocs = pop_readlocs(s.chanlocs); 
% chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
% chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
% load(s.chanlocs_lw);

cmap = rs_prepareTopo;

%% select files 
subs = {'SubM', 'SubT'};
temp = listdlg('ListString', subs, 'SelectionMode', 'single');
sub = subs{temp}; clear temp 

files = dir; files = {files.name}; % files( ~endsWith(files, '.set')) = [];
files( ~endsWith(files, [ sub '.set'])) = []; % select subject
if isempty(files), error 'RS error: no files found'; end
sel = listdlg('ListString', files, 'SelectionMode','multiple', 'ListSize',[400,600], 'InitialValue',1:length(files));
% sel = 2; warning( '! BYPASSING SESSION SELECTION' )
files = files(sel);
nfiles = length(files);

%% loop through sessions
seshlist = [];
for f = 1 :nfiles
    
    % get subject ID and session
    temp = strsplit(files{f}(1:end-4)); sesh = temp{end-1}(1:3); sub = temp{end}; clear temp

    %% load file 
    EEG = pop_loadset(files{f}); 
    
    % select only EEG channels
    EEG = pop_select(EEG, 'channel', find( strcmp({EEG.chanlocs.type},'EEG') & ~startsWith({EEG.chanlocs.labels},'EXG')) );

    % store for later merge
    EEG_all(f) = EEG;
    

end % end session loop

% merge files
EEG_all = pop_mergeset(EEG_all, 1:length(EEG_all));
disp 'FINISHED MERGING FILES'


%% filter data stuff for use in loop below

%    % low-pass filter lfp to remove line-noise and other frequency signals, using same lowpass as EEG
% %     LFP = pop_eegfiltnew( LFP, 'hicutoff', EEG.etc.preproc.freqs(2)); % low-pass
%     LFP_all = pop_eegfiltnew( LFP_all, 'hicutoff', 40); % low-pass (higher band of 40 just in case freqs are higher)


    %% condition loop 
    disp 'running PLV ...'

    % get info 
    nconds = length(unique({EEG_all.event.type}));
    for cond = 1:nconds
    if any(strcmp({EEG_all.event.type},s.conds{cond})) % check this cond actually exists for this session first
        % split by condition
        EEG = pop_select(EEG_all, 'trial', find(strcmp({EEG_all.event.type},s.conds{cond})) );
        ntrials = size(EEG.data,3); 

        % rename EEG variable
        EEG.setname = [s.conds{cond} ' mergedSesh ' sub];
            
        % run PLV analysis
        cfg = [];
        cfg.bl_win = [ -1 0 ];
        edges = nonLinspace(1,30 , 5, 'exp10' );
        for k = 1:length(edges)-1
            cfg.filt_bands(k,:) = [ edges(k), edges(k+1) ];
        end
        PLV = rs_eeg_phaseLockingValue(EEG, cfg); 

        % export to lw  
        for f = 1:length(PLV)
            cd([ s.savePath filesep 'lw'])
            rs_convert_lab2lw_V1( PLV(f)  , PLV(f).setname, chanlocs_lw );
            fprintf('... saved file: %s\n', [ cd filesep PLV(f).setname ])
        end

        
        fprintf('\ncompleted cond %d/%d\n', cond, length(s.conds))
    end % end if statement checking whether this cond exists for this session
    end % end cond loop
    fprintf('\n~~~ FINISHED ALL CONDS')
    


%% END
% cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

