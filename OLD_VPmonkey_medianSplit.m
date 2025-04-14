% 
% 
% 
%       cluster MUA/LFP responses & plot histograms of their depth values
% 
% 
%         !! todo:  CURRENTLY this fails with visual because peaks are funny ... very rough 
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

file_types   = {'MUA','EEG'};
split_type   = 'EEG';
split_chan   = 'Cz'; 
split_range  = [0, 0.3]; % range to search for peaks within
split_peakwin = 5/1000; % width either side around each peak to average around in seconds
split_thresh = 3; % min amplitude to find a peak 
split_npeaks = 3;

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
        if strcmp(file_types{ft},'EEG')

            files = dir; files = {files.name}; 
            files = files( endsWith(files,[ sub '.set']) & ...
                ~contains(files,'all') )';
            files = files(  contains(files,'merged_EEG') );

            % check number of files match the session list
            if length(files2load) ~= nfiles
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
            data.(file_type)(f) = pop_loadset('filename',files2load{f},'filepath',s.savePath);
        end
        fprintf('%s ... finished loading\n', file_type)
        
    %% END loop through file-types2
    end

    %% extract data channel to be used for splitting
    split_chan_ind = find(strcmpi(  {data.(split_type)(1).chanlocs.labels}, split_chan   ));
    for f = 1:nfiles
        data.split(f) = pop_select( data.(split_type)(f), 'channel', split_chan_ind);
    end

    %% low-pass filter data used to compute median splits
    for f = 1:nfiles
        data.split(f) = pop_eegfiltnew( data.split(f), 'hicutoff', 30 );
    end

    %% get only certain time-range of data
    for f = 1:nfiles
        data.split(f) = pop_select( data.split(f), 'time',  split_range);
    end
    

    %% LOOP THROUGH CONDITIONS
    for c = 1:length(s.conds)

        %% concatenate data to be used for splitting
        d = [];
        for f = 1:nfiles
            temp = data.split(f).data;
            ents = {data.(split_type)(f).event.type};
            temp(:,:,~strcmp(ents,s.conds{c})) = []; % remove epochs which correspond to different modality
            temp = permute(temp,[3 2 1]);
            % ? could do z-scoring within session here
            d = [d; temp];
        end

        %% use average of all trials to find peak latencies
        dm = mean(d);
        [~,locs] = findpeaks(abs(dm),"MinPeakHeight",split_thresh, 'NPeaks',split_npeaks)
        peakwin = round(split_peakwin * data.split(1).srate);

        %% loop through peaks and use to median split all data
        for pk = 1:split_npeaks

            loc = locs(pk);

            data4split = mean( d(:, loc-peakwin:loc+peakwin) ,2); 
            %       figure; histogram( data4split )

            % median split
            [~, inds] = sort( data4split );
            mid = floor(length( inds)/2);

            %% loop through file types and median split
            for ft = 1:length(file_types)
                file_type = file_types{ft};

                % append all trials of this file type
                merged = pop_mergeset( data.(file_type), 1:nfiles );

                % get this condition
                merged = pop_select(merged,'trial', find(strcmp({merged.event.type},s.conds{c})));

                % split 
                for m = 1:2 % low to high
                    if m == 1
                        dtemp  = pop_select(merged,'trial', inds(1:mid) ); 
                        msplitcode = 'M_LOW';
                    else
                        dtemp  = pop_select(merged,'trial', inds(mid+1:end) ); 
                        msplitcode = 'M_HIGH';
                    end

                    % export
                    % downsample for efficiency
                    if ~strcmp( file_type, 'MISC' )
                        dtemp  = pop_rs_downsample(dtemp,s.dsf);
                    end
        
                    % if LFP or MUA, just select average for simplicity
                    if strcmp(file_types{ft},'MUA') 
                        dtemp = pop_select(dtemp, 'channel', find(strcmp({dtemp.chanlocs.labels},'MUA_avg')) );
                    end
                    if strcmp(file_types{ft},'LFP')
                        dtemp = pop_select(dtemp, 'channel', find(strcmp({dtemp.chanlocs.labels},'LFP_avg')) );
                    end
                    
                    % crop for lw
                    dtemp = pop_select(dtemp, 'time', s.xlims.lw.(file_type) );
        
                    % replace all NaNs with zeros for easier lw plotting
                    dtemp.data(isnan(dtemp.data)) = 0;
            
                    % export to lw
                    cd([s.savePath filesep 'lw'])
                    saveName = [ 'ep_' s.conds{c} ...
                        ' medsplit_' split_type ' pk_' num2str(pk)  ...
                        ' ' msplitcode ...
                        ' merged_' file_type ' ' sub ];
                    rs_convert_lab2lw_V1( dtemp ...
                        , saveName, [] );

                end % median split loop
               


            %% END LOOP THROUGH FILE TYPES
            end 

        %% END LOOP THROUGH PEAKS
        end 
        
    %% END LOOP THROUGH CONDITIONS
    end


%% END loop through subs
end

%%
cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

