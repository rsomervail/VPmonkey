% 
% 
%       
%       compute hilbert envelope of gamma power
% 
% 
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
% homedir = ['/media/rick/Rick_LabData_3/neuro/iannettilab_ongoing/VPmonkey/data/clean_' subfold '_OLD'];
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
s.savePath = homedir;


% plot limits 
% xlims
s.xlims.plot.EEG  = [-0.2 0.6];  
s.xlims.plot.LFP  = [-0.2 0.6];  
s.xlims.plot.MUA  = [-0.2 0.6];  
s.xlims.plot.EYE  = [-1   6]; 
s.xlims.plot.MISC = [-0.2 0.5];
% ylims - amplitudes
s.ylims.plot.EEG  = [-20 20];  
s.ylims.plot.LFP  = [-80 80];   
s.ylims.plot.MUA  = [-2 2];  
s.ylims.plot.EYE  = [-2e-2   2e-2]; 
% ylims - tvals
s.ylims.tvals.EEG  = [-40 40];  
s.ylims.tvals.LFP  = [-80 80];  
s.ylims.tvals.MUA  = [-40 40]; 
s.ylims.tvals.EYE  = [-40 40];
% topoplot clims - amplitudes
s.clims.SubM.AUD = [-5 5];
s.clims.SubM.SOM = [-5 5];
s.clims.SubM.VIS = [-5 5];
s.clims.SubT.AUD = [-5 5];
s.clims.SubT.SOM = [-5 5];
s.clims.SubT.VIS = [-5 5];

% EEG topo peaks  ? this will find the closest peak on the average to the specified points
s.topopeaks.SubM.AUD = [30 95 160]/1000;  %!!!!!!!!!! REDO WHEN HAVE FINAL AVERAGES
s.topopeaks.SubM.SOM = [36 100 180] / 1000
s.topopeaks.SubM.VIS = [24 45  67  94 120 290 ]/1000;
s.topopeaks.SubT.AUD = [30 95 160]/1000;
s.topopeaks.SubT.SOM = [36 100 180] / 1000
s.topopeaks.SubT.VIS = [24 45  67  94 120 290 ]/1000;

s.dsf = 2; % downsample before lw export

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

s.conds = {'AUD','SOM','VIS'};

str_exclude = {}; % exclude all file names containing these strings

str_include = {}; % include only file names containing these strings


subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};


s.bandpass = [50 100];
% s.bandpass = [50 300];


%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    savestr = []; % tracking which modules are used to reject trials

    %% get list of files
    cd(homedir)
    files_all = dir; files_all = {files_all.name}; 

    % filter files
    files_all = files_all(  endsWith(files_all,[ sub '.set']) & ...
        contains(files_all,'LFP') & ~contains(files_all,'all') )';
    if ~isempty(str_exclude)
        files_all( cellfun( @(x)  any( contains(x,str_exclude)) , files_all)  ) = [];
    end
    if ~isempty(str_include)
        files_all( cellfun( @(x)  any(~contains(x,str_include)) , files_all)  ) = [];
    end
    files = files_all;
    nfiles = length(files);

    %% loop through files
    for f = 1:nfiles
        %% load LFP
        file2load = files{f};
        LFP = pop_loadset('filename',file2load,'filepath',homedir, 'loadmode',1:5);
        
    
        %% band-pass filter
        LFP = pop_eegfiltnew( LFP, 'locutoff', s.bandpass(1) , 'hicutoff', s.bandpass(2));
    
        %% make dummy EEGLAB structure
        GAM = LFP; 
        GAM.data = nan(size(GAM.data));
        for c = 1:GAM.nbchan
            GAM.chanlocs(c).labels = strrep(GAM.chanlocs(c).labels,'LFP','GAM');
        end
    
        %% compute Hilbert envelope
        %
        for c = 1:LFP.nbchan
            for e = 1:LFP.trials
                temp = double(LFP.data(c,:,e));
                GAM.data(c,:,e) = abs(hilbert(temp));  % ? use second parameter window length?  round(LFP.srate * 20/1000)
            end
        end
        
        %% baseline correction
        GAM = pop_rmbase(GAM,[-400 -100]);
    
        %% compute mean of channels
        GAM.data(end+1,:,:) = mean(GAM.data,'omitnan');
        GAM.chanlocs(end+1) = GAM.chanlocs(end);
        GAM.chanlocs(end).labels = 'GAM_avg';
        GAM.nbchan = GAM.nbchan + 1;
    
        %% save in EEGLAB
        filename = strrep(file2load,'merged_LFP','merged_GAM');
        filename = [ 'BP_' num2str(s.bandpass(1)) '_' num2str(s.bandpass(2)) 'Hz ' filename];
        pop_saveset(GAM,'filepath', homedir, 'filename',filename);

    %% END loop through files
    end
    
%% END loop through subs
end 

%%
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'
