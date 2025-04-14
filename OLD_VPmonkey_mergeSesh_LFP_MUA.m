% 
% 
% 
%       splits electrodes per session and saves them individually so they can be considered together
%       across sessions by depth
%       
%           
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

s.saveas = 'allElec'; % name for this set of sessions SHOULD ALWAYS INCLUDE "all" at the start for later string matching


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

% file_types = {'LFP', 'MUA', 'GAM'};
% file_types = {'LFP', 'MUA', 'GAM', 'EEG'};
file_types = {'MUA', 'EEG'};


%% import intracranial information table
file_depths = [getRoot '/VPmonkey/SUA_protocol.xlsx'];
tbl_depths = VPmonkey_import_electrodeDepths( file_depths  );

%% get sessions (assumes the necessary files exist)
sesh = tbl_depths.sesh;
nfiles = length(sesh); % ? should not filter sessions at this stage .. can do that more easily later when analysing

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
            files = files(  contains(files,'yesICA') ); % CLEANED WITH ICA
%             files = files( ~contains(files,'yesICA') ); % NOT cleaned with ICA

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
                data.(file_type)(f) = pop_loadset('filename',files2load{f},'filepath',homedir ...
                    , 'loadmode', 1:5); % only load actual electrodes, not avg of electrodes
            else
                data.(file_type)(f) = pop_loadset('filename',files2load{f},'filepath',homedir);
            end
        end
        fprintf('%s ... finished loading\n', file_type)
        
    %% END loop through file-types
    end

    %% get depths of each loaded electrode
    % get row indices of selected sessions
    inds = cellfun( @(x) find(strcmp(x, tbl_depths.sesh))  , sesh );

    % get common field names across electrodes
    flds = {'Dura','TOA','Depth','DepthRelTOA','DepthRelDura'};

    % loop through sessions and stack depths of each electrode in order
    tbl_depths_byElec = [];
    for f = 1:nfiles
        
        % loop through electrodes
        for k = 1:5

            % loop through fields
            tbl_temp = table;
            for fd = 1:length(flds)
                temp = tbl_depths.([flds{fd} '_' sub(end) num2str(k)]);
                tbl_temp.(flds{fd}) = temp(inds(f));
            end
            tbl_depths_byElec = [tbl_depths_byElec; tbl_temp];
        end
    end
    save( [s.savePath filesep 'tbl_depths_byElec_' sub],'tbl_depths_byElec')


    %% LOOP THROUGH FILE TYPES AND SPLIT BY ELECTRODE
    for ft = 1:length(file_types)
        file_type = file_types{ft};

        %% LOOP THROUGH FILES
        count = 1;
        for f = 1:nfiles

            %% if not intracranial then obviously just repeat the same data 5 times
            if strcmp(file_type,'EEG') % ? could do this also with EYE but not sure what the point would be

                % make 5 copies of the same EEG data because the channels match
                for k = 1:5   
                    dtemp = data.(file_type)(f);
                    dtemp.filename = [ 'elec_' sprintf('%03d',count) ' ' dtemp.filename];
                    dtemp.setname = dtemp.filename(1:end-4);
                    pop_saveset(dtemp, 'filename', dtemp.filename, 'filepath',s.savePath);

                    count = count + 1; % increment electrode counter
                end
    
            %% if NOT EEG then split electrodes
            else

                % loop through electrodes and save data separately
                for k = 1:5
                    dtemp = pop_select(data.(file_type)(f),'channel', k);
                    dtemp.filename = [ 'elec_' sprintf('%03d',count) ' ' dtemp.filename];
                    dtemp.setname = dtemp.filename(1:end-4);
                    pop_saveset(dtemp, 'filename', dtemp.filename, 'filepath',s.savePath);

                    count = count + 1; % increment electrode counter
                end


            %% END IF EEG OR NOT
            end

        %% END FILE LOOP
        end

    %% END file type loop
    end 

    
    


%% END loop through subs
end

%%
cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

