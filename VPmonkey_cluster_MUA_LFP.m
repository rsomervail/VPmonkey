% 
% 
% 
%       cluster MUA/LFP responses & plot histograms of their depth values
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

s.saveas = 'allSesh'; % name for this set of sessions SHOULD ALWAYS INCLUDE "all" at the start for later string matching

% epoching 
s.xlims.lw.EEG  = [-0.5 1];  
s.xlims.lw.LFP  = [-0.5 1];  
s.xlims.lw.MUA  = [-0.5 1];  
s.xlims.lw.GAM  = [-0.5 1];  
s.xlims.lw.MISC = [-0.2 0.5];
s.xlims.lw.EYE  = [-0.5 6]; 

s.dsf = 2; % downsample before lw export

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

file_types = {'GAM'};

clust_type = 'GAM';

% number of clusters 
nclust = 2;

% binmetric = 'Depth';
% binmetric = 'DepthRelDura';
binmetric = 'DepthRelTOA';

%% import intracranial information table
file_depths = [getRoot '/VPmonkey/SUA_protocol.xlsx'];
tbl_depths = VPmonkey_import_electrodeDepths( file_depths  );

%% select sessions (assumes the necessary files exist)
sesh = tbl_depths.sesh;
% sel = listdlg('ListString',sesh,'SelectionMode','multiple');
% sesh = sesh(sel);
nfiles = length(sesh);

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};
    clear ntrials_all

    %% loop through file types & load selected data 
    for ft = 1:length(file_types)
        file_type = file_types{ft};
    
        %% get session filenames
        cd(homedir)
        if strcmp(file_types{ft},'EEG') || strcmp(file_types{ft},'GAM')

            files = dir; files = {files.name}; 
            files = files( endsWith(files,[ sub '.set']) & ...
                ~contains(files,'all') & ~contains(files,'bin') )';
            files = files(  contains(files,['merged_' file_types{ft}]) );

            % check number of files match the session list
            if length(files) ~= nfiles
                error 'MISMATCH IN FILE-NUMBERS, MISSING SESSIONS OR NEED TO HANDLE FILE-FILTERING?';
            end

            % check that these are the correct sessions
            temp = cell2mat(regexp( files, 'S\d\d' ));
            if length(unique(temp)) ~= 1, error 'files with same format should have session code in the same place...'; end
            temp2 = cellfun( @(x) x(temp(1):temp(1)+2) ,files ,'Uniformoutput',false);
            if ~isequal(temp2,sesh)
                error 'file session codes don''t match session list'
            end

            files2load = files;

        else
            files2load = repmat({['merged_' file_type ' S##_MT ' sub '.set']}, length(sesh), 1);
            files2load = cellfun( @(x,y) strrep(x,'S##',y) ,files2load, sesh, 'UniformOutput', false);
        end
        
        %% load sessions for this file_type
        for f = 1:nfiles
            if strcmp(file_type,'MUA') || strcmp(file_type,'LFP') || strcmp(file_type,'GAM')
                data.(file_type)(f) = pop_loadset('filename',files2load{f},'filepath',s.savePath ...
                    , 'loadmode', 1:5); % only load actual electrodes, not avg of electrodes
            else
                data.(file_type)(f) = pop_loadset('filename',files2load{f},'filepath',s.savePath);
            end
        end
        fprintf('%s ... finished loading\n', file_type)
        
    %% END loop through file-types2
    end

    %% LOOP THROUGH CONDITIONS
    for c = 1:length(s.conds)

        %% stack all channels to be used for clustering
        d = [];
        for f = 1:nfiles
            temp = data.(clust_type)(f).data;
            ents = {data.(clust_type)(f).event.type};
            temp(:,:,~strcmp(ents,s.conds{c})) = []; % remove epochs which correspond to different modality
            temp = mean(temp,3);
            d = [d; temp];
        end


        %% get depths of each loaded electrode
        inds = cellfun( @(x) find(strcmp(x, tbl_depths.sesh))  , sesh );
        deps = [];
        for k = 1:5
            temp = tbl_depths.([binmetric '_' sub(end) num2str(k)]);
            deps = [ deps, temp(inds) ];
        end
        deps_r = reshape(deps',[size(deps,1) * size(deps,2),1]);
            
    
        %% remove zero or nan channels
        inds2remove = any(isnan(d),2);
        inds2remove = any(d==0,2) | inds2remove;
        d(inds2remove,:) = [];
        deps_r(inds2remove) = [];
    
        %% prepare data for clustering
        % crop
        times = data.(clust_type)(1).times/1000;
        t1 = findnearest(times,0);
        t2 = findnearest(times,0.3);
    
        d2 = d(:,t1:t2);
        times2 = times(t1:t2);
    
        % z-score 
        d2 = (d2 - mean(d2,2)) ./ std(d2,[],2);
    
        %% cluster
    %     [ids] = kmeans( d2, nclust, 'Distance', 'correlation');
        [ids] = kmeans( d2, nclust );
    
        counts = histcounts(ids);
    
        % plot cluster waveforms & t-tests 
        for k = 1:nclust
            dtemp = d(ids==k,:);
            [~,pvals] = ttest(dtemp);
    
            figure('name', [sub ' ' s.conds{c} ' _ cluster ' num2str(k) ] ,'NumberTitle','off');
            subplot(2,1,1); 
            plot(times, mean(dtemp) ); xlim([-0.2 0.6])
            title(['cluster ' num2str(k) ': n = ' num2str(counts(k))])
    
            subplot(2,1,2); 
            plot(times, pvals ); xlim([-0.2 0.6])
        end
    
        % plot histograms of electrode depths for each cluster 
        figure('name', [sub ' ' s.conds{c}] ,'NumberTitle','off');
        for k = 1:nclust
            deps_r_clust{k} = deps_r(ids==k);
            histogram( deps_r_clust{k}, 15, ...
                'DisplayName',['cluster ' num2str(k) ]); hold on
            legend
        end

    %% END LOOP THROUGH CONDITIONS
    end


%% END loop through subs
end

%%
cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

