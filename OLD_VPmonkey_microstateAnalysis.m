%
%    run microstate analysis on each condition separately using mixed gaussian models & compute stats (e.g. frequency of state transitions etc)
% 
% 
%%
clc
clearvars
close all

subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
% subfold = 'main_example';

homedir = ['/media/rick/Rick_LabData/neuro/iannettilab_ongoing/VPmonkey/data/clean_' subfold];
% homedir = ['/home/rick/neuro/iannettilab/Projects/VP monkey/data/clean_' subfold];
cd(homedir);

addpath('/home/rick/neuro/iannettilab/Projects/VP monkey/Giac_ToolBox')
addpath('/home/rick/neuro/iannettilab/Projects/VP monkey/scripts')

% try 
%     eeglab
% catch
%     addlab
% end

%% SETTINGS
s = [];
s.savePath =  ['/media/rick/Rick_LabData/neuro/iannettilab_ongoing/VPmonkey/data/clean_' subfold]; % mkdir(s.savePath)

s.conds = {'VIS', 'SOM', 'AUD'};
s.epwin = [ -2  2 ]; % window of interest


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
    EEG_all = pop_loadset(files{f}); 
    
    % select only EEG channels
    EEG_all = pop_select(EEG_all, 'channel', find( strcmp({EEG_all.chanlocs.type},'EEG') & ~startsWith({EEG_all.chanlocs.labels},'EXG')) );

    % get info
    if f == 1
        nchans  = EEG_all.nbchan;
        chanlocs = EEG_all.chanlocs;
    end

%    % low-pass filter lfp to remove line-noise and other frequency signals, using same lowpass as EEG
% %     LFP = pop_eegfiltnew( LFP, 'hicutoff', EEG.etc.preproc.freqs(2)); % low-pass
%     LFP_all = pop_eegfiltnew( LFP_all, 'hicutoff', 40); % low-pass (higher band of 40 just in case freqs are higher)

    % crop epochs to focus on ERP
    EEG_all = pop_select(EEG_all, 'time', s.epwin);
    
    % get info 
    nconds = length(unique({EEG_all.event.type}));

    %% condition loop 
    disp 'running microstate analysis ...'
    for cond = 1:nconds
    if any(strcmp({EEG_all.event.type},s.conds{cond})) % check this cond actually exists for this session first
        % split by condition
        EEG = pop_select(EEG_all, 'trial', find(strcmp({EEG_all.event.type},s.conds{cond})) );
        ntrials = size(EEG.data,3); 
            
        % run microstate analysis
        cfg = [];
        cfg.nstates = 3:10;
        cfg.regval  = 0.2;
        cfg.std_thresh = 10;
        EEG = rs_eeg_microstateAnalysis(EEG, cfg); 

        
        fprintf('\ncompleted cond %d/%d\n', cond, length(s.conds))
    end % end if statement checking whether this cond exists for this session
    end % end cond loop
    fprintf('\n~~~ completed file %d/%d\n\n\n',f,length(files))
    

end % end session loop
disp 'FINISHED FILES'


%% END
% cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

