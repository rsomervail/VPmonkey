% 
% 
% 
%       cross-correlate intra-cranial data of different modalities (e.g. MUA x LFP)
%           (uses split electrode datasets (byElectrode)
% 
%       ? correlation with MUA gives nothing, probs because of poor sampling/high single trial noise
%       
%       ? 1 drawback here is that it's only convenient to compare diff signals from the same electrode
%           ... trial-by-trial correlations only make sense within session
%           ... but can maybe still find an interesting effect of depth?
%           ... or at the very least just find the overall relationship between MUA, LFP and GAM
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

pathlw  = [getRoot 'toolbox' filesep 'letswave7-master'];
pathlab = [getRoot 'toolbox' filesep 'eeglab2022.1'];

%% switch to eeglab
evalc 'addpath(genpath(pathlab));';
evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''fieldtrip-20210311'' ]))'; % remove fieldtrip to tidy things up a bit
evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''Biosig3.7.9'' ]))'; % remove biosig to tidy things up a bit
evalc 'rmpath(genpath(pathlw));';

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
% file_type_comps = {'LFP', 'MUA'}; 
% file_type_comps = {'LFP', 'GAM'}; 
% file_type_comps = {'EEG', 'MUA'}; % EEG MUST BE FIRST  
% file_type_comps = {'EEG', 'LFP'}; 
file_type_comps = {'EEG', 'GAM'}; 
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
            temp = pop_loadset('filename',files2load{f},'filepath', s.savePath);

            % if EEG then run lowpass filter
            if strcmp(file_type,'EEG')
                temp = pop_eegfiltnew(temp, 'hicutoff', 30 );
            end

            % downsample data 
            temp = pop_rs_downsample(temp, 8);
            
            data.(file_type)(f) = temp;

        end
        fprintf('%s ... finished loading\n', file_type)
        
    %% END loop through file-types
    end


    %% LOOP THROUGH COMPARISONS
    tIN = tic;
    for ft = 1:size(file_type_comps,1)

        rvals_all = cell(nfiles,length(s.conds));
        pvals_all = cell(nfiles,length(s.conds));

        
        %% LOOP THROUGH ELECTRODE FILES
        for f = 1:nfiles

            %% get the two data files (one session, one electrode, from each modality)
            d1 = data.(file_type_comps{ft,1})(f);
            d2 = data.(file_type_comps{ft,2})(f);

            %% check comparison pair have same epochs
            if d1.trials ~= d2.trials || ~isequal({d1.event.type},{d2.event.type})
                error('TRIAL NUMBER MISMATCH \nd1: %s\nd2: %s',d1.filename,d2.filename);
            end
            enttypes = {d1.event.type};

            %% LOOP THROUGH CONDITIONS
            for cond = 1:length(s.conds)

                %% get trials from this condition
                inds = find(strcmp(enttypes,s.conds{cond}));
                d1cond = pop_select(d1,'trial',inds);
                d2cond = pop_select(d2,'trial',inds);

                %% check whether one file is EEG
                if ~any(strcmp(file_types,'EEG'))
                    dat1 = squeeze(d1cond.data)';
                    dat2 = squeeze(d2cond.data)';

                    %% CROSS-CORRELATION

                    [rvals,pvals] = corr(dat1,dat2, 'rows','complete');
    
                    figure; imagesc(rvals); 
                    xlabel(file_type_comps{ft,2}); ylabel(file_type_comps{ft,1});
                    ax = gca;
                    ax.YDir = 'normal';
                    % !!! Mike Cohen book has some nice code for replacing the ticks here
    %                 ax.XTick = 
                    
                    figure; imagesc(pvals,[0,0.05]); 
                    xlabel(file_type_comps{ft,2}); ylabel(file_type_comps{ft,1});
                    ax = gca;
                    ax.YDir = 'normal';


                else % EEG
                    dat1 = squeeze(d1cond.data); % assuming EEG is dat1
                    dat2 = squeeze(d2cond.data)';

                    nchans = size(dat1,1); % assuming EEG is dat1
                    npnts_eeg = size(dat1,2); 
                    npnts_lfp = size(dat2,2); 
                    rvals = single(nan(nchans, npnts_eeg, npnts_lfp));
                    pvals = rvals;

                    %% CROSS-CORRELATION
                    % loop through channels
                    parfor c = 1:size(dat1,1)

                        
                        [rvals(c,:,:),pvals(c,:,:)] = ...
                            corr( squeeze(dat1(c,:,:))',dat2, 'rows','complete');
    

%                     figure; imagesc(rvals); 
%                     xlabel(file_type_comps{ft,2}); ylabel(file_type_comps{ft,1});
%                     ax = gca;
%                     ax.YDir = 'normal';
%                     % !!! Mike Cohen book has some nice code for replacing the ticks here
%     %                 ax.XTick = 
%                     
%                     figure; imagesc(pvals,[0,0.05]); 
%                     xlabel(file_type_comps{ft,2}); ylabel(file_type_comps{ft,1});
%                     ax = gca;
%                     ax.YDir = 'normal';


                    end % chan loop

                    rvals_all{f,cond} = rvals;
                    pvals_all{f,cond} = pvals;


                end

            %% output
            fprintf('finished file %03d/%03d cond %d/%d in %.2f mins\n\n' ...
                        , f, nfiles, cond, length(s.conds), toc(tIN)/60 )
      
            %% END CONDITION LOOP
            end
            

        %% END ELECTRODE FILE LOOP
        end

        mkdir([cd '/lw']); cd([cd '/lw'])

        %% switch to letswave
        evalc 'rmpath(genpath(pathlab));';
        evalc 'addpath(genpath(pathlw));';

        %% LOOP THROUGH CONDITIONS
        for cond = 1:length(s.conds)

            %% REFORMAT FOR LW EXPORT
            rvals = cat(4, rvals_all{:,cond} );
%             pvals = cat(4, pvals_all{:,cond} );

            %% REJECT NAN EPOCHS
            badsesh = find(any(isnan(rvals),[1,2,3]));
            rvals(:,:,:,badsesh) = [];
%             pvals(:,:,:,badsesh) = [];

            %% EXPORT TO LW
            
            % export to lw
            if any(strcmp(file_types,'EEG'))
                option.dimension_descriptors={'channels','X','Y','epochs'};
            else
                error 'HAVENT CODED EXPORT FOR NO EEG YET'
            end
            option.unit='amplitude';
            option.xunit='time';
            option.yunit='time';
            option.xstart   = d1.times(1)/1000;
            option.xstep    = 1/d1.srate;
            option.ystart   = d2.times(1)/1000;
            option.ystep    = 1/d2.srate;
            option.is_save  = 1;
            
            filename = [   'crossCorr rvals goodSesh ' s.conds{cond} ' ' ...
                file_type_comps{ft,1} '_' file_type_comps{ft,2} ' ' sub  ];
            option.filename = filename;
            FLW_import_mat.get_lwdata(rvals,option);

            % rename EEG channels
            if any(strcmp(file_types,'EEG'))
                option = struct('filename',filename);
                lwdata = FLW_load.get_lwdata(option);
                option = struct( ...
                    'old_channel',{{'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','C21','C22','C23','C24','C25','C26','C27','C28'}},...
                    'new_channel',{{'O2','Oz','O1','PO4','POz','PO3','P4','P2','P1','P3','CP4','CP2','CPz','CP1','CP3','CZ','FC6','FC4','FC2','FCz','FC1','FC3','FC5','F2','Fz','F1','AF4','AF3'}}, ...
                    'suffix','','is_save',1);
                lwdata = FLW_electrode_labels.get_lwdata(lwdata,option);
            end

            % compute average
            fprintf('computing average ...\n')
            option = struct('filename',filename);
            lwdata = FLW_load.get_lwdata(option);
            option = struct('operation','average','suffix','avg','is_save',1);
            lwdata = FLW_average_epochs.get_lwdata(lwdata,option);

            % compute t-ttest
            fprintf('computing t-test ...\n')
            option = struct('filename',filename);
            lwdata = FLW_load.get_lwdata(option);
            option = struct('constant',0,'tails','both','alpha',0.05,'suffix','ttest','is_save',1);
            lwdata = FLW_ttest_constant.get_lwdata(lwdata,option);

            %
            cond

        %% END COND LOOP
        end

    %% END COMPARISON LOOP
    end

%% END loop through subs
end

%%
cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

