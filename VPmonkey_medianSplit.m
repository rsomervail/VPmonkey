% 
% 
% 
%       median split of data of N types according to values of another type
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

pathlw  = [getRoot 'toolbox' filesep 'letswave7-master'];
pathlab = [getRoot 'toolbox' filesep 'eeglab2022.1'];

%% SWITCH TO EEGLAB
evalc 'addpath(genpath(pathlab));';
evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''fieldtrip-20210311'' ]))'; % remove fieldtrip to tidy things up a bit
evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''Biosig3.7.9'' ]))'; % remove biosig to tidy things up a bit
evalc 'rmpath(genpath(pathlw));';

%% SETTINGS
s = [];
s.savePath =  [homedir filesep 'lw']; % mkdir(s.savePath)
s.savePath_figs =  [ getRoot '/VPmonkey/paper/figures/raw' ]; % mkdir(s.savePath_figs)

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

% file_types   = {'LFP','EYE','MUA','EEG'};
file_types   = {'EYE','LFP'};
% split_type   = 'EEG'; split_chan   = 'Cz'; 
split_type   = 'LFP'; split_chan = 'LFP'; % ? makes the most sense analytically but haven't tested script out like this yet
split_range  = [0, 0.3]; % range to search for peaks within
split_peakwin = 4/1000; % width either side of each single-trial peak to average around in seconds
split_peaksearchwin = 25/1000; % width either side of avg peak to search for single-trial peak within 
split_thresh = 0.1; % min z-score to find a peak  ? has to be weirdly small atm, probs because z-scoring is done across all conds
split_npeaks = 3;

%% import intracranial information table
tbl_depths = VPmonkey_import_electrodeDepths;

%% loop through subjects
for sb = 1:length(subs)
    sub = subs{sb};

    %% SWITCH TO EEGLAB
    evalc 'addpath(genpath(pathlab));';
    evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''fieldtrip-20210311'' ]))'; % remove fieldtrip to tidy things up a bit
    evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''Biosig3.7.9'' ]))'; % remove biosig to tidy things up a bit
    evalc 'rmpath(genpath(pathlw));';

    %% LOAD SELECTED DATA
    % by session (EEG)
    cfg = struct;
    cfg.sub = sub;
        % trial-rejection criteria
        ar = struct;
        ar.method = 'time';
        ar.metric = 'median';
        ar.thresh = 3;
        ar.timeprop = 0.15; % ?
        ar.chanprop = 0.15; % ?
        cfg.autoar = ar;
    cfg.filetypes = file_types;
    cfg.zscore    = file_types;
    cfg.zscore_cond = true;
    cfg.byElectrode = true;
    cfg.average   = false;
    cfg.mergesesh = true;
    cfg.filt_preproc = 'allSesh'; % 'allGoodEEG', 'allPerfEEG', 'allSesh'
    cfg.EEG_include = {'LP_30 ged_horiz','yesICA'};
    cfg.EEG_exclude = {'noICA'};
    data = VPmonkey_mergeSesh(cfg);
    
    fprintf('... finished loading data\n')

    %% extract data channel to be used for splitting
    split_chan_ind = find(strcmpi(  {data.(split_type).chanlocs.labels}, split_chan   ));
    data.split = pop_select( data.(split_type), 'channel', split_chan_ind);

    %% low-pass filter data used to compute median splits
    data.split = pop_eegfiltnew( data.split, 'hicutoff', 30 );

    %% high-pass filter data used to compute median splits
    data.split = pop_eegfiltnew( data.split, 'locutoff', 1 );

    %% copy filtered version to it's main field for later export of appropriately filtered signal
    data.(split_type) = data.split;

    %% get only certain time-range of data
    data.split = pop_select( data.split, 'time',  split_range);

    %% LOOP THROUGH CONDITIONS
    data.split = pop_rs_removeNonZeroEvents(data.split);
    for c = 1:length(s.conds)

        %% concatenate data to be used for splitting
        d = data.split.data;
        ents = {data.split.event.type};
        d(:,:,~strcmp(ents,s.conds{c})) = []; % remove epochs which correspond to different modality
        d = permute(d,[3 2 1]); % trials, samples, channels (? nchans=1)

        %% use average of all trials to find peak latencies
        % find peaks
        dm = mean(d);
        [pkstemp,locs] = findpeaks(abs(dm),"MinPeakHeight",split_thresh)
        peakwin = round(split_peakwin * data.split.srate);

        % find biggest N peaks
        [pkstemp, inds] = sort( pkstemp,'descend')
        locs = locs(inds(1:split_npeaks));
        locs = sort(locs);
        
        % get peak signs
        pksigns = sign(dm(locs));

        % plot peaks on average
        figure; 
        plot(dm); hold on;
        for pk = 1:split_npeaks
            plot([locs(pk),locs(pk)],ylim,'k')
        end
        title([sub ' - peaks 4 split - ' s.conds{c}  ])
        % save fig to matlab format
        figname = [sub '_medsplit_' split_type '_peaks4split_' s.conds{c}  ];
        temp = gcf;
        savefig(temp, [s.savePath_figs filesep figname])
        % export as .eps file
        temp.Renderer = 'painters';
        exportgraphics(temp, [s.savePath_figs filesep figname '.eps'])

        % compute peak amplitudes for each trial and each peak (NEW - trial-by-trial peak detection)
        pks = nan(size(d,1),split_npeaks);
        searchwin = round((split_peaksearchwin*data.split.srate));
        for pk = 1:split_npeaks
            loc = locs(pk);

            % loop through trials
            for trl = 1:size(d,1)

                % get trial
                dtrial = d(trl,:);
%                 figure; plot(dtrial)

                % get search window around mean peak latency
                dsearchwin = dtrial(loc-searchwin:loc+searchwin);
%                 figure; plot(dsearchwin)

                % if this is a negative peak, flip the sign
                if pksigns(pk)==-1
                    dsearchwin = -dsearchwin;
                end

                % find peak in single-trial search window
                [amps,loc2] = findpeaks(dsearchwin);

                % if a single-trial peak is found within this range, then find its location within the whole trial
                if ~isempty(loc2)

                    % if multiple peaks returned then sort them and return the largest peak
                    if length(loc2) > 1
                        [~,indstemp] = sort(amps,'descend');
                        loc2 = loc2(indstemp(1));
                    end

                    % convert this index to it's equivalent in the whole trial
                    loc3 = loc2 + (loc-searchwin-1);

                else % else if no single-trial peak is found within range, then just use the avg peak for this trial
                    loc3 = loc;
                    warning('no single-trial peak found for cond %s pk %d trl %d',s.conds{c},pk,trl)
                end

                % extract mean signal around single-trial peak
                pks(trl,pk) = mean( dtrial(loc3-peakwin:loc3+peakwin)); 

%                 figure; plot(dtrial); hold on; 
%                 plot([loc, loc],ylim,'g')
%                 plot([loc3,loc3],ylim,'b')
%                 plot([loc3-peakwin,loc3-peakwin],ylim,'k')
%                 plot([loc3+peakwin,loc3+peakwin],ylim,'k')
                
            end % trial loop

            %       figure; histogram( pks(:,pk) )
        end

%         % OLD - fixed peak latency for each trial
%         compute peak amplitudes for each trial and each peak
%         pks = nan(size(d,1),split_npeaks);
%         for pk = 1:split_npeaks
%             loc = locs(pk);
%             pks(:,pk) = mean( d(:, loc-peakwin:loc+peakwin) ,2); 
%             %       figure; histogram( pks(:,pk) )
%         end

        if any(isnan(pks),'all')
            error 'NAN PEAKS FOUND??'
        end

        % NEW - add sum of sign-flipped peaks as another "virtual peak"
        pktemp = pks .* pksigns;
        pks = [pks,  sum(pktemp,2) ];
        pksigns = [pksigns, 1];
        %       figure; histogram( pks(:,end) )

%         % OLD - add SD across peaks as another "virtual peak"
%         pks = [pks, std(pks,[],2)];
%         pksigns = [pksigns, 1];
%         %       figure; histogram( pks(:,end) )

        %% loop through peaks and use to median split all data
        for pk = 1:split_npeaks+1

            % median split
            if pksigns(pk) == 1 % positive peak
                [~, inds] = sort( pks(:,pk) );
            elseif pksigns(pk) == -1 % negative peak
                [temp, inds] = sort( pks(:,pk), 'descend' );
            end
            mid = floor(length(inds)/2);

            %% loop through file types and median split
            for ft = 1:length(file_types)
                file_type = file_types{ft};

                % SWITCH TO EEGLAB
                evalc 'addpath(genpath(pathlab));';
                evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''fieldtrip-20210311'' ]))'; % remove fieldtrip to tidy things up a bit
                evalc 'rmpath(genpath([ pathlab filesep ''plugins'' filesep ''Biosig3.7.9'' ]))'; % remove biosig to tidy things up a bit
                evalc 'rmpath(genpath(pathlw));';

                % get this condition
                dtemp = data.(file_types{ft});
                dtemp = pop_rs_removeNonZeroEvents(dtemp);
                dtemp = pop_select(dtemp,...
                    'trial', find(strcmp({dtemp.event.type},s.conds{c})));

                % split 
                for m = 1:2 % low to high
                    if m == 1
                        dtemp_bin  = pop_select(dtemp,'trial', inds(1:mid) ); 
                        msplitcode = 'M_LOW';
                    else
                        dtemp_bin  = pop_select(dtemp,'trial', inds(mid+1:end) ); 
                        msplitcode = 'M_HIGH';
                    end

                    % EXPORT TO LW
                    % downsample for efficiency
                    if ~strcmp( file_type, 'MISC' )
                        dtemp_bin  = pop_rs_downsample(dtemp_bin,s.dsf);
                    end
        
                    % crop for lw
                    dtemp_bin = pop_select(dtemp_bin, 'time', s.xlims.lw.(file_type) );
        
                    % replace all NaNs with zeros for easier lw plotting
                    dtemp_bin.data(isnan(dtemp_bin.data)) = 0;
            
                    % export to lw
                    if pk==(split_npeaks+1)
                        pkstr = 'pk_SUM';
                    else
                        pkstr = ['pk_' num2str(pk)];
                    end
                    cd(s.savePath)
                    saveName = [ 'ep_' s.conds{c} ...
                        ' medsplit_s' split_type '_' strjoin(file_types(~strcmp(file_types,split_type)),'_')  ...
                        ' ' pkstr  ...
                        ' merged_' file_type ...
                        ' ' msplitcode ...
                         ' ' sub ];
                    if strcmp(file_type,'EEG')
                        rs_convert_lab2lw_V1( dtemp_bin ...
                            , saveName, chanlocs_lw);
                    else
                        rs_convert_lab2lw_V1( dtemp_bin ...
                            , saveName, []);
                    end

                end % median split loop

                %% SWITCH TO LETSWAVE
                evalc 'rmpath(genpath(pathlab));';
                evalc 'addpath(genpath(pathlw));';

                %% LW - COMPUTE AVERAGE AND T-TEST
                % average - LOW
                saveName_LOW = strrep(saveName,'M_HIGH','M_LOW');
                option_load = struct('filename',[saveName_LOW '.lw6']);
                lwdata = FLW_load.get_lwdata(option_load);
                option = struct('operation','average','suffix','avg','is_save',1);
                lwdata = FLW_average_epochs.get_lwdata(lwdata,option);

                % average - HIGH
                option_load = struct('filename',[saveName '.lw6']);
                lwdata = FLW_load.get_lwdata(option_load);
                option = struct('operation','average','suffix','avg','is_save',1);
                lwdata = FLW_average_epochs.get_lwdata(lwdata,option);
            
                % t-test across subjects (cluster perm)
                option = struct('filename',{{[saveName_LOW '.lw6'], [saveName '.lw6']}});
                lwdata = FLW_load.get_lwdataset(option);
                option = struct('test_type','two-sample','tails','both',...
                    'ref_dataset',1, 'show_progress',1,'suffix','ttest','is_save',1);
                lwdata = FLW_ttest.get_lwdataset(lwdata,option);
%                 % t-test across subjects (cluster perm)
%                 option = struct('filename',{{[saveName_LOW '.lw6'], [saveName '.lw6']}});
%                 lwdata = FLW_load.get_lwdataset(option);
%                 option = struct('test_type','two-sample','tails','both',...
%                     'ref_dataset',1,'alpha',0.05,'permutation',1,'cluster_threshold',0.05,...
%                     'num_permutations',2000,'show_progress',1,'suffix','ttest','is_save',1);
%                 lwdata = FLW_ttest.get_lwdataset(lwdata,option);

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

