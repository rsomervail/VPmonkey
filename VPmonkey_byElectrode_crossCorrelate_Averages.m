% 
% 
% 
%       cross-correlate intra-cranial data of different modalities (e.g. MUA x LFP)
%           (uses split electrode datasets (byElectrode)
%       
%       ? 1 drawback here is that it's only convenient to compare diff signals from the same electrode
%           ... trial-by-trial correlations only make sense within session
%           ... but can maybe still find an interesting effect of depth?
%           ... or at the very least just find the overall relationship between MUA, LFP and GAM
%
%       ! probably need to do this across trials and also across electrodes/sessions using averages
% 
%%
clc
clearvars
close all

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

saveas = '';
% saveas = '_SPK_thresh_45';

homedir = ['/media/rick/Rick_LabData_3/neuro/iannettilab_ongoing/VPmonkey/data/clean_' subfold saveas];
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
s.savePath =  [homedir '/byElectrode/']; if ~exist(s.savePath,'dir'), mkdir(s.savePath); end

% s.saveas = '';

% channel locations for these datasets
s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
s.chanlay  = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
chanlocs = pop_readlocs(s.chanlocs); 
chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
load(s.chanlocs_lw);

s.conds = {'AUD','SOM', 'VIS'};

subs = {'SubM','SubT'};
% subs = {'SubM'};
% subs = {'SubT'};

% choose pairs of data types to correlate with each other
file_type_comps = {'LFP', 'MUA'}; 
% file_type_comps = {'LFP', 'GAM'}; 
% file_type_comps = {'LFP', 'MUA'; 'MUA', 'GAM'; 'LFP', 'GAM'; }; % pairs of data types to correlate with each other

file_types = unique(file_type_comps); % no need to change this manually

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    clear ntrials_all

    %% load tbl with depths per electrode
    load( [s.savePath filesep 'tbl_depths_byElec_' sub])
    
    %% loop through file types & load selected data 
    for ft = 1:length(file_types)
        file_type = file_types{ft};
    
        %% get filenames
        cd(s.savePath)
    
        files = dir; files = {files.name}; 
        files = files( ...
             endsWith(files,[ sub '.set']) ...
            & startsWith(files,'elec_')  ... 
            & contains(files,file_type) ... 
            )';
        

%         files() = []; % ? could have filters here for different preproc settings

        files2load = files;
        nfiles = length(files);
        
        %% load sessions for this file_type
        for f = 1:nfiles
            data.(file_type)(f) = pop_loadset('filename',files2load{f},'filepath', s.savePath);
        end
        fprintf('%s ... finished loading\n', file_type)
        
    %% END loop through file-types
    end

    % !! once have data nicely formatted above, can loop through file_types vs file_types 
%               and do channel-by-channel correlations across time

    %% LOOP THROUGH COMPARISONS
    for ft = 1:size(file_type_comps,1)

        %% LOOP THROUGH CONDITIONS
        for cond = 1:length(s.conds)
        
            %% LOOP THROUGH ELECTRODE FILES
            dat1 = nan( data.(file_type_comps{ft,1})(1).pnts ,nfiles);
            dat2 = nan( data.(file_type_comps{ft,2})(1).pnts, nfiles);
            for f = 1:nfiles
    
                %% get the two data files (one session, one electrode, from each modality)
                d1 = data.(file_type_comps{ft,1})(f);
                d2 = data.(file_type_comps{ft,2})(f);
    
                %% check comparison pair have same epochs
                if d1.trials ~= d2.trials || ~isequal({d1.event.type},{d2.event.type})
                    error('TRIAL NUMBER MISMATCH \nd1: %s\nd2: %s',d1.filename,d2.filename);
                end
                enttypes = {d1.event.type};

                %% get trials from this condition
                inds = find(strcmp(enttypes,s.conds{cond}));
                d1cond = pop_select(d1,'trial',inds);
                d2cond = pop_select(d2,'trial',inds);

                %% compute averages across trials
                temp1 = squeeze(d1cond.data);
                temp2 = squeeze(d2cond.data);
                
                temp1 = mean(temp1,2);
                temp2 = mean(temp2,2);

                % store averages for later
                dat1(:,f) = temp1;
                dat2(:,f) = temp2;

                

            %% END ELECTRODE LOOP
            end

            %% reformat data
            dat1 = dat1';
            dat2 = dat2';

            %% compute grand averages
            dat1m = mean(dat1,'omitnan');  
            dat2m = mean(dat2,'omitnan'); 

            times1 = data.(file_type_comps{ft,1})(1).times/1000;
            times2 = data.(file_type_comps{ft,2})(1).times/1000;

            figure; 
            plot(times1,dat1m)
            yyaxis right
            plot(times2,dat2m)

            %% correlate 
            [rvals, pvals] = corr(dat1,dat2, 'rows', 'complete'); 

            %% plot results % !! do this properly and add saving and exporting

            figure; imagesc(rvals); 
            xlabel(file_type_comps{ft,2}); ylabel(file_type_comps{ft,1});
            ax = gca;
            ax.YDir = 'normal';
            % !!! Mike Cohen book has some nice code for replacing the ticks here
%                 ax.XTick = 

            alpha = 0.05;
            rvals_thresh = rvals; rvals_thresh(pvals < alpha) = 0;
            figure; imagesc(rvals_thresh); 
            xlabel(file_type_comps{ft,2}); ylabel(file_type_comps{ft,1});
            ax = gca;
            ax.YDir = 'normal';
            % !!! Mike Cohen book has some nice code for replacing the ticks here
%                 ax.XTick = 
            
            figure; imagesc(pvals,[0,0.05]); 
            xlabel(file_type_comps{ft,2}); ylabel(file_type_comps{ft,1});
            ax = gca;
            ax.YDir = 'normal';

        %% END ELECTRODE FILE LOOP
        end

    %% END COMPARISON LOOP
    end

%% END loop through subs
end

%%
cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

