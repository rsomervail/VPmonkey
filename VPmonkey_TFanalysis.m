%
%
%       run wavelet time-frequency analysis on monkey VP data (EEG & LFP)
% 
%           - use merged EEG if doing EEG, otherwise use byElectrode for LFP 
%           - compute the average per session and then merge these across sessions/electrodes etc
% 
% 
%           !! MIGHT NOT HAVE ENOUGH BASELINE ... NEED TO RE-RUN PREPROC SPECIFICALLY FOR TF??
%
%           - Richard Somervail 2023
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

pathlw  = [getRoot 'toolbox' filesep 'letswave7-master'];
pathlab = [getRoot 'toolbox' filesep 'eeglab2022.1'];

% try 
%     eeglab
% catch
%     addlab
% end

%% SETTINGS
s = [];
s.savePath =  homedir;

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

s.conds = {'AUD','SOM','VIS'}; 
nconds = length(s.conds);

subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};

file_types = {'EEG'};
% file_types = {'EEG','LFP'};

s.eeg_filter = 'allGoodEEG';
% s.lfp_filter = '';

%% LOAD TABLE WITH EEG PREPROC INFO 
tbl_preproc = VPmonkey_load_table_preproc;

%% LOOP SUBS
for sb = 1:length(subs)
    sub = subs{sb};

%% 
cd(homedir)
files_all = dir; files_all = {files_all.name}'; 
files_all = files_all(   endsWith(files_all, [ sub '.set']) );
files_all = files_all(  ~contains(files_all,'all') );

tin = tic;

%% loop through filetypes
for ft = 1:length(file_types)

    %% switch to eeglab
    evalc 'addpath(genpath(pathlab));';
    evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''fieldtrip-20210311'' ]))'; % remove fieldtrip to tidy things up a bit
    evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''Biosig3.7.9'' ]))'; % remove biosig to tidy things up a bit
    evalc 'rmpath(genpath(pathlw));';


    %% get files
    switch file_types{ft}
        case 'EEG'
    
            % filter files
            files2load = files_all(  endsWith(files_all,[ sub '.set']) & ...
                startsWith(files_all,'yesICA')  );


            % FILTER BY EEG DATA QUALITY

            % filter by subject
            tbl = tbl_preproc;
            tbl = tbl( strcmp(tbl.sub,sub(end)) ,:);
            if length(files2load) ~= size(tbl,1), error 'FILE NUMBER MISMATCH!!!'; end
        
            % filter by bad EEG sessions
            switch s.eeg_filter
                case 'allGoodEEG' % good in at least 2 conditions
                    temp = sum( [ tbl.rej_aud tbl.rej_som tbl.rej_vis ], 2 );
                    indy = temp <= 1; % i.e. no more than 1 bad condition
                    goodsesh = tbl.sesh(indy);
        
                case 'allPerfEEG' % good in all 3 conditions
                    temp = sum( [ tbl.rej_aud tbl.rej_som tbl.rej_vis ], 2 );
                    indy = temp == 0; % i.e. no bad conditions
                    goodsesh = tbl.sesh(indy);
        
                case 'allSesh'
                    goodsesh = tbl.sesh;
            end

            % filter by good sessions
            indy = cellfun( @(x) contains(x,goodsesh), files2load );
            files2load = files2load(indy);
            nfiles = length(files2load);
    
        case 'LFP'
    
    
    end

%% LOOP THROUGH FILES
for f = 1:nfiles

    %% load file
    EEG = pop_loadset(files2load{f});

    %% LOOP conditions
    for cond = 1:nconds

        %% get this condition
        EEG_cond = pop_select( EEG, 'trial', find(strcmpi({EEG.event.type},s.conds{cond})) );

        %% RUN WAVELET ANALYSIS   ! time to code my own one and make it really efficient, use parallel processing
    
    %     numcycles = 6;
        numcycles = [3 10]; % range between lower and upper frequencies analysed
    
        power = cell(EEG.nbchan,1);
        times = cell(EEG.nbchan,1);
        freqs = cell(EEG.nbchan,1);
        parfor c = 1:EEG.nbchan
            [~,~,~,times{c},freqs{c},~,~,tfdata] ...
                = pop_newtimef(EEG_cond, 1, c, [] ...
                    ,numcycles, 'freqs',[5 100] ...
                    ,  'nfreqs', 40 ...
                    ,  'freqscale', 'linear' ... %  log   linear
                    , 'timesout', EEG.pnts/4 ... % ? this will need to change if data format changes
                    , 'baseline', nan ...
                    , 'plotersp','off', 'plotitc','off'...
                    , 'verbose', 'off'); 

            power{c,1} = abs(tfdata);  
            fprintf('.')
        end
        fprintf('\n')
        power = cat(4,power{:});
        freqs = freqs{1};
        times = times{1};
        nfreqs = length(freqs); 
        nsamps = length(times);
         
        power = permute(power, [2,1,3,4]);

        %% BASELINE CORRECT
        blind = [ findnearest(-120, times) , findnearest(-50,times)   ];
        blind = blind(1):blind(2);

        % step 1 - single trial z-scoring
        power_bl = (power - mean(power)) ./ std(power); 

        % step 2 - baseline z-score the average across trials
        temp = squeeze(mean(power_bl,3));
        power_bl_mean = (temp - mean(temp(blind,:,:))) ./ std(temp(blind,:,:)); 

        %% STORE MEAN FOR MERGING ACROSS SESSIONS
        if ~exist('power_bl_mean_all','var')
            power_bl_mean_all = nan( [ nfiles nconds size(power_bl_mean) ]  );
        end
        power_bl_mean_all(f,cond,:,:,:) = power_bl_mean;

    %% END CONDITION LOOP
    fprintf('cond %d/%d completed - %s\n',cond,nconds,files2load{f})
    end

%% END FILE LOOP  %!!! ADD TIME 
fprintf('file %d/%d complete in %.2f mins - %s\n\n',f,nfiles, toc(tin)/60, files2load{f})
end

%% EXPORT TO LW
% switch to letswave
evalc 'rmpath(genpath(pathlab));';
evalc 'addpath(genpath(pathlw));';

cd([homedir '/lw'])
for cond = 1:nconds
 
    dtemp = squeeze(power_bl_mean_all(:,cond,:,:,:));

    %% import data (TF)
    option.dimension_descriptors={'epochs','X','Y','channels'};
    option.unit     =   'amplitude';
    option.xunit    =   'time';
    option.yunit    =   'frequency';
    option.xstart   = times(1)/1000;
    option.xstep    = mean(diff(times))/1000;
    option.ystart   = freqs(1);
    option.ystep    = mean(diff(freqs));
    option.is_save  = 1;
    
    cycstr = ['cyc_' strrep(num2str(numcycles),' ','_')];
    option.filename = ['TF bl_gc ' cycstr ' ' s.conds{cond} ' merged_' file_types{ft} ' ' sub ];
    FLW_import_mat.get_lwdata( dtemp,option);

    filename = option.filename;

    %% ADD CHANNEL LABELS
    option = struct('filename',filename);
    lwdata = FLW_load.get_lwdata(option);
    option = struct('old_channel',{{'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28'}},'new_channel',{{'O2','Oz','O1','PO4','POz','PO3','P4','P2','P1','P3','CP4','CP2','CPz','CP1','CP3','CZ','FC6','FC4','FC2','FCz','FC1','FC3','FC5','F2','Fz','F1','AF4','AF3'}},'suffix','','is_save',1);
    lwdata = FLW_electrode_labels.get_lwdata(lwdata,option);

    %% ADD CHANNEL LOCATIONS  
    option = struct('filename',filename);
    lwdata= FLW_load.get_lwdata(option);
    option=struct('filepath','/home/rick/Dropbox/Somervail_Richard/toolbox/letswave7-master/res/electrodes/spherical_locations/monkeyEEG_LW.locs','suffix','','is_save',1);
    lwdata= FLW_electrode_location_assign.get_lwdata(lwdata,option);

    %% AVERAGE across sessions
    option = struct('filename',filename);
    lwdata = FLW_load.get_lwdata(option);
    option = struct('operation','average','suffix','avg','is_save',1);
    lwdata = FLW_average_epochs.get_lwdata(lwdata,option);

    %% T-TEST ACROSS SESSIONS
    option = struct('filename',filename);
    lwdata = FLW_load.get_lwdata(option);
    option = struct('constant',0,'tails','both','alpha',0.05,'suffix','ttest','is_save',1);
    lwdata = FLW_ttest_constant.get_lwdata(lwdata,option);

end % cond loop




%% END LOOP THROUGH FILETYPES
end


%% END SUB LOOP
end

