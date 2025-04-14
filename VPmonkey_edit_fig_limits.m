

%% SELECT AND LOAD FIGURE
clearvars
close all
figpath = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/paper/figures/raw';
cd(figpath)

files = dir;
files = files(~[files.isdir]);
files = {files.name}';
files = files(endsWith(files,'.fig'));
files = cellfun(@(x) strrep(x,'.fig',''), files, 'UniformOutput', false);

sel = listdlg('ListString',files,'ListSize',[600,800], 'SelectionMode','single');
file = files{sel};

fig = open([file '.fig']);

%% X-LIMITS
% xlim([   ]);
figure(fig);

%% Y-LIMITS
% ylim([ -40  40  ]); % LFP t-values
ylim([ -2  2  ]); % pupil averages
figure(fig);

%% EXPORT FIGURE
figure(fig);
rs_saveExportFig(fig, figpath, file);



