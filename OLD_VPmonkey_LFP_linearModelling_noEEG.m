% 
% 
%   ** unfinished **
% 
%       predict LFP amplitude with simple intercept LME models + including depth and recording site
% 
%       ? doing with matlab LME for now but eventually do student t in R
% 
%%
clc
clearvars
close all


homedir = '/media/rick/Rick_LabData_3/neuro/datasets/VPmonkey/data/clean_main';
cd(homedir);

addpath([getRoot '/VPmonkey/Giac_ToolBox'])
addpath([getRoot '/VPmonkey/scripts'])

figsdir = '/home/rick/Dropbox/Somervail_Richard/VPmonkey/paper/figures';

subs = {'SubM','SubT'};

filetype = 'LFP';

%%
tic
for sb = 1:length(subs)
    sub = subs{sb};

    % artifact rejection settings
    ar = struct;
    ar.method = 'time';
    ar.metric = 'median';
    ar.thresh = 3;
    ar.timeprop = 0.1; % ? too strict?
    ar.chanprop = 0.1; % ? too strict?

    % load data
    cfg = struct; cfg.exportlw = false;
    cfg.sub = sub;
    cfg.filetypes = {'LFP'};
    cfg.LFP_include = {'LP_30'};
%     cfg.LFP_exclude = {'LP_30'};
    cfg.byElectrode = true;
    cfg.average   = false;
    cfg.mergesesh = false;
    cfg.zscore = {}; 
%     cfg.zscore = cfg.filetypes; 
%     cfg.zscore_win = [-0.2 0.6]; 
%     cfg.zscore_cond  = true;
    [data,tbl_depths] = VPmonkey_mergeSesh(cfg);
    LFP = data.LFP;

    sites = unique(tbl_depths.dist_angle);
    nsites = length(sites);
    nplots = find( ((1:6).^2) >= nsites ,1,'first'); % for subplots later

    % 
    nelec = length(LFP);
    for el = 1:nelec
        LFP(el) = pop_rs_downsample(LFP(el),4);
        LFP(el) = pop_select(LFP(el),'time',[-0.2,0.6]);
    end

    %% EXTRACT DATA
    
    % EXTRACT ALL TRIALS 
    ents = cell(nelec,1); 
    d  = cell(nelec,1);
    for el = 1:nelec
        tempents = LFP(el).event;
        for e = 1:length(tempents)
            tempents(e).elec = el;
        end
        ents{el} = tempents;
        d{el} = squeeze(LFP(el).data);
    end
    ents = cat(2,ents{:});
    d = cat(2,d{:});
    times = LFP(1).times/1000;
    nsamps = length(times);

    conds  = unique({ents.type});
    nconds = length(conds);

    %% RUN MODELS
    for cond = 1:nconds
    
        %% extract this condition
        inds = strcmp({ents.type},conds{cond});
        ents_cond = ents(inds);
        dcond = d(:,inds);

        %% LOOP THROUGH TIMEPOINTS
%         dpred = nan(size(dcond));
        site  = arrayfun(@(x) tbl_depths.dist_angle(x), [ents_cond.elec])';
        Depth = arrayfun(@(x) tbl_depths.Depth(x), [ents_cond.elec])';
        Depth = (Depth-mean(Depth))./std(Depth);
        tbl = table(site,Depth);
        tbl.site = categorical(tbl.site);

        for t = 1:nsamps
            
            %% format data for LME matrix function
            y = double(dcond(t,:)');
            tbl.y = y;
%             x = []; %!! HOW TO FORMAT

            %
%             m = fitlmematrix()


            % TEMP INEFFICIENT
            models = {'y ~ 1',...
                'y ~ 1 + Depth',...
                'y ~ -1 + (1 | site)'...
                'y ~ -1 + Depth + (1 | site)'...
                'y ~ -1 + (Depth | site)'};
            nmodels = length(models);
            m = cell(nmodels,1);
            BIC = nan(nmodels,1);
%             for mi = 1:nmodels
                mi = 2;% !!!
                m{mi} = fitlme(tbl, models{mi}, 'DummyVarCoding', 'full' );
                BIC(mi) = m{mi}.ModelCriterion.BIC;
%             end
%             figure; bar(BIC); ax = gca; ax.XTickLabels = models;

            est(t) = m{2}.Coefficients.Estimate(2);
            ci(t) = m{2}.Coefficients.Upper(2) - m{2}.Coefficients.Estimate(2);

%             m = fitlm(tbl,'y ~ -1 + site','DummyVarCoding', 'full');
            % ! can contrast this with the LMEs above using BIC

%             dpred(t,:) = m.predict; % ? only really interesting when split by site
    
        % !!  from LM can get CI for site-wise averages

            t
        end


    %% END CONDITION LOOP
    end
%%
end % sub loop 
