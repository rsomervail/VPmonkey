% 
% 
%       Plot LFP electrode averages colour-coded by depth
%           - Richard Somervail, 2025
%
%   !! UNFINISHED
% 
%%
clc
clearvars
close all

% subfold = 'main'; % main experiment data
% subfold = 'pilots_techtests';
subfold = 'main';

homedir = ['/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_' subfold];
cd(homedir);

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

sub = 'SubM';
% sub = 'SubT';

%% by electrode (LFP, MUA etc)
cfg = struct; cfg.exportlw = false;
% cfg.autoar = ar;
cfg.sub = sub;
cfg.filetypes = {'LFP'};
cfg.LFP_include = {'LP_30'};
cfg.byElectrode = true;
cfg.average   = true;
cfg.mergesesh = true;
cfg.zscore = {}; 
% cfg.zscore_win = [-0.2 0.6];
% cfg.zscore_cond = true;
[temp, tbl_depths] = VPmonkey_mergeSesh(cfg);
LFP = temp.LFP;

%% plot LFP electrode averages colour-coded by electrode depth

cond2plot = 'AUD';
% cond2plot = 'SOM';
% cond2plot = 'VIS';

% extract data
inds = find(strcmpi({LFP.event.type},cond2plot));
LFP2plot = pop_select(LFP,'trial',inds);
tbl_depths_cond = tbl_depths(inds,:);
d = squeeze(LFP2plot.data);
times = LFP2plot.times/1000;

% % plot individual electrodes
% figure; plot(times,d); xlim([-0.2 0.4]); 

% plot mean across electrodes
dm = mean(d,2); figure; plot(times,dm); xlim([-0.2 0.4])

%% plot scatterplot within window against electrode depth
depthmetric = 'Depth';
% depthmetric = 'DepthRelTOA';
% depthmetric = 'DepthRelDura';
% indy = findnearest(times,0.020):findnearest(times,0.055); % N WAVE
indy = findnearest(times,0.072):findnearest(times,0.116); % P WAVE
d2plot = mean(d(indy,:));
x = tbl_depths_cond.(depthmetric);
y = d2plot';
y(isnan(x)) = []; x(isnan(x)) = [];
figure; scatter( x , y ); xlabel(depthmetric); ylabel 'LFP (uV)'; lsline
[r,p] = corr(x,y)

%% plot r across time
% !! note this is not the same as slope which would be better probably 
%    (can probably plot coefficient timecourse w/ CI)

depthmetric = 'Depth';
% depthmetric = 'DepthRelTOA';
% depthmetric = 'DepthRelDura';

x = tbl_depths.(depthmetric);
y = d';

[r,p] = corr(x,y);

xlims = [-0.2 0.4];

figure; 
nrows = 2;
subplot(nrows,1,1); plot(times,mean(y)); xlim(xlims);
subplot(nrows,1,2); plot(times,r); xlim(xlims);





