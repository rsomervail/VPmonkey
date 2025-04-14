% 
% 
%   ** unfinished **
% 
%       plot LFP/MUA etc by recording site
% 
%       ? eventually should arrange these geometrically and superimpose over a MK brain map top down view
% 
%%
clc
clearvars
close all


homedir = '/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_main';
cd(homedir);

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

figsdir = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/paper/figures/raw';

subs = {'SubM','SubT'};

filetype = 'LFP';

%%
tic
for sb = 1:length(subs)
    sub = subs{sb};

    % by electrode (LFP, MUA etc)
    cfg = struct; cfg.exportlw = false;
    cfg.sub = sub;
    cfg.filetypes = {filetype};
    cfg.LFP_include = {'LP_30'};
%     cfg.LFP_exclude = {'LP_30'};
    cfg.byElectrode = true;
    cfg.average   = true;
    cfg.mergesesh = true;
    cfg.zscore = {};  % ? z-scoring may obscure important differences in VP presence!
%     cfg.zscore_win = [-0.2 0.6];
%     cfg.zscore_cond  = true;
    [data,tbl_depths] = VPmonkey_mergeSesh(cfg);

    % extract data for each condition
    ents = data.LFP.event;
    conds_all = {ents.type};
    conds = unique(conds_all);
    nconds = length(conds);
    clear LFP
    for cond = 1:nconds
        inds = find(strcmp(conds_all,conds{cond}));
        LFP(cond) = pop_select(data.LFP,'trial',inds);
        tbl_depths_cond{cond} = tbl_depths(inds,:);
    end
    
    sites = unique(tbl_depths.dist_angle);
    nsites = length(sites);
    nplots = find( ((1:6).^2) >= nsites ,1,'first'); % for subplots later


    times = LFP(1).times/1000;

    %% PLOT
    for cond = 1:nconds

        sites_all = tbl_depths_cond{cond}.dist_angle;
        d = squeeze(LFP(cond).data);
    
        fig_title = ['plot-bySite ' filetype ' '  sub ' ' conds{cond}];
        fig = figure('name',fig_title,'numbertitle','off');
        for k = 1:nsites
            subplot(nplots,nplots,k);
            
            inds = find(strcmp(sites_all, sites{k}));
            if ~isempty(inds)
    
                dtemp = d(:,inds);
                
                plot(times,dtemp); hold on;
                plot(times,mean(dtemp,2),'k','LineWidth',2)
                
                xlim([-0.2 0.6])
                ylim([-100,100])
%                 ylim([-3,3])

            end

            title(strrep(sites{k},'_',' '))

        end % site loop

        % EXPORT
        rs_saveExportFig(fig, figsdir, fig_title,1);

    end % condition loop

end % sub loop 
