%
%    merge LFPs by depth bin and export to lw
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
s.savePath =  homedir;

s.conds = {'VIS', 'SOM', 'AUD'};
s.epwin = [ -0.2  0.6 ]; % window of interest

% % channel locations for these datasets
% s.chanlocs = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW.locs';
% s.chanlay  = '/home/rick/neuro/iannettilab/Projects/VP monkey/Giac_ToolBox/Layouts/Layout_Monkey_EEG.mat';  % fieldtrip version of chanlocs for Giac interpolation
% s.chanlocs_lw = '/home/rick/Dropbox/Somervail_Richard/EEG Configurations/monkeyEEG_LW_chanlocs.mat'; % for lw export
% chanlocs = pop_readlocs(s.chanlocs); 
% chanlocs = chanlocs(end:-1:1);  % reordering because order is reversed in locs file...
% chanlocs(16).labels = 'CZ'; % also renaming Cz to CZ as in biosemi
% load(s.chanlocs_lw);

cmap = rs_prepareTopo;

file_types = {'MUA','LFP'};
% file_types = {'MUA','LFP','EEG'};

%% import intracranial information table
file_depths = [getRoot '/VPmonkey/SUA_protocol.xlsx'];
tbl_depths = VPmonkey_import_electrodeDepths( file_depths  );

%% select subject 
subs = {'SubM', 'SubT'};
temp = listdlg('ListString', subs, 'SelectionMode', 'single');
sub = subs{temp}; clear temp 

%% select files 
files = dir; files = {files.name}; 
files( ~endsWith(files, [ sub '.set'])) = []; % select subject
files( ~endsWith(files, [ sub '.set'])) = []; % select subject
if isempty(files), error 'RS error: no files found'; end
sel = listdlg('ListString', files, 'SelectionMode','multiple', 'ListSize',[400,600], 'InitialValue',1:length(files));
% sel = 2; warning( '! BYPASSING SESSION SELECTION' )
files = files(sel);
nfiles = length(files);

%% loop through sessions
seshlist = [];
tbl = [];
for f = 1 :nfiles
    
    % get subject ID and session
    temp = strsplit(files{f}(1:end-4)); sesh = temp{end-1}(1:3); sub = temp{end}; clear temp

    %% load LFP 
    if ~exist( 'LFP_inds', 'var')
        LFP = pop_loadset(files{f}, cd ); 
        LFP_inds = find( strcmp({LFP.chanlocs.type},'LFP'));
        LFP = pop_select(LFP, 'channel', LFP_inds );
        nlfp = length(LFP_inds);
    else
        LFP = pop_loadset(files{f}, cd, LFP_inds ); 
    end

%     % low-pass filter lfp to remove line-noise and other frequency signals, using same lowpass as EEG
% %     LFP = pop_eegfiltnew( LFP, 'hicutoff', EEG.etc.preproc.freqs(2)); % low-pass
% %     LFP_all = pop_eegfiltnew( LFP_all, 'hicutoff', 40); % low-pass (higher band of 40 just in case freqs are higher)
%     % !! could also do time-frequency envelope/power of e.g. 55 - 200 Hz gamma here

    % z-score to normalize ! needs to be done manually :/
%     LFP = pop_rmbase(LFP, );
    
    % loop through LFP electrodes in this file
        fprintf([ repmat('.',1,nlfp) '\n\n' ] )
        for l = 1:nlfp

            % check if LFP has depth value, indicating that it was not broken
            if isnan( tbl_depths.(['Depth_' sub(end) num2str(l)])(strcmp(tbl_depths.sesh, sesh)) )  % check if LFP electrode was broken, if so then skip
                continue   
            end

            % extract only this one channel
            LFP_singleChan = pop_select(LFP, 'channel', l);
            LFP_singleChan.chanlocs(1).labels = 'LFP';
            
            if ~exist('LFP_all','var')
                LFP_all = LFP_singleChan;
            else
                LFP_all(end+1) = LFP_singleChan; 
            end

            tbl(end+1  ,1).Dura = tbl_depths.(['Dura_' sub(end) num2str(l)])(strcmp(tbl_depths.sesh, sesh)); 
            tbl(end    ,1).Depth = tbl_depths.(['Depth_' sub(end) num2str(l)])(strcmp(tbl_depths.sesh, sesh)); 
            tbl(end    ,1).TOA = tbl_depths.(['TOA_' sub(end) num2str(l)])(strcmp(tbl_depths.sesh, sesh)); 
            tbl(end    ,1).DepthRelTOA  = tbl(end  ,1).Depth - tbl(end  ,1).TOA;
            tbl(end    ,1).DepthRelDura = tbl(end  ,1).Depth - tbl(end  ,1).Dura;

        
            fprintf('\b|\n')
        end % end LFP loop
            
        indy_lfp = ((f-1)*nlfp + 1) : ((f-1)*nlfp + nlfp) ; % indices for this set of 5 LFPs

    fprintf('\n~~~ completed file %d/%d\n\n\n',f,length(files))
    
    seshlist = [seshlist; repmat({sesh}, nlfp,1)]; 
end % end session loop
disp 'FINISHED EXPORTING DATA'

%% export files to lw but labelled by bin

% make bins
depth_measure = 'DepthRelDura';
% depth_measure = 'DepthRelTOA';
nbins = 3;
[ bin_counts, bin_edges,  bindy ] = histcounts(  [tbl.(depth_measure)]  , nbins );

% loop through LFP electrodes and export each one according to 
for l = 1:length(LFP_all)

    if bindy(l) > 0
        LFP = LFP_all(l);
    
        % export to lw
        cd([ s.savePath filesep 'lw'])
        saveName = [ 'bin_'  num2str(bindy(l)) ' LFP_' num2str(l) ' binMeas_' depth_measure ' onlyLFP '  sub  ];
        rs_convert_lab2lw_V1( LFP  , saveName, LFP.chanlocs );
        fprintf('... saved LFP_%d: %s\n', l, [s.savePath filesep saveName]) 
    end

end


%% END
% cd(s.savePath)
disp '~~~~~~~~~~FINISHED~~~~~~~~~~'

