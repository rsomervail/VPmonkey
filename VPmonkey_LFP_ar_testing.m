%
%
%       testing the ar_lfp trial rejection in VPmonkey_mergeSesh
% 
%       1) thresh_abs = [200,5/100; 400,0]
%           - loses no good trials, but misses a few mid amplitude, high variance bad trials (e.g. SubT 20 & 65)
%           - also misses the odd brief spike that comes close to 400 uV
% 
%       2) thresh_abs = [150,10/100; 300,0]
%           - still misses some of those high variance, mid amplitude trials, but improves
%             things a lot without substantially rejecting good trials (e.g. SubT f20 is pretty good, f65 not great) 
%             SubT 121 has one bad trial still
% 
%       3) thresh_abs = [150,5/100; 300,0]
%           - a little aggressive, catches a handful of big VPs with peaks at 200 uV in early SubM files
%           - also somehow misses some bad trials in e.g. SubM f70-74
%           - but cleans up SubT f20 very nicely! f21 and f65 less so..
%           - if I really need very low noise rate then maybe this is ok, but othrewise option 2 is better
% 
%      ** best all-rounder so far is option 2 **
%       
%      cfg.ar_lfp.maxbad 
%           - 0.25 only results in one electrode being rejected (from SubT, none in SubM)
%           - 0.10 results in 4 being rejected (again all from SubT)
%           
%       ** 0.25 is a good default, 0.10 would be stricted and probably remove the really bad electrodes
%          with trials not caught by the absolute threshold criterion above **
% 
% 
%%
for sb = 1:2

    if sb == 1
        sub = 'SubM';
    else
        sub = 'SubT';
    end

    % by electrode (LFP, MUA etc)
    cfg = struct; cfg.exportlw = false;
    cfg.sub = sub;
%     cfg.ar_lfp.thresh_abs = [200,5/100; 400,0]; % first column is absolute thresholds, second is max proportion of timepoints above threshold
    cfg.ar_lfp.thresh_abs = [150,10/100; 300,0]; % first column is absolute thresholds, second is max proportion of timepoints above threshold
%     cfg.ar_lfp.thresh_abs = [150,5/100; 300,0]; % first column is absolute thresholds, second is max proportion of timepoints above threshold
    cfg.ar_lfp.maxbad = 0.25; % max number of epochs to reject before rejecting electrode
    cfg.filetypes = {'LFP'};
    cfg.LFP_include = {'LP_30'};
%     cfg.LFP_exclude = {'LP_30'};
    cfg.byElectrode = true;
    cfg.average   = false;
    cfg.mergesesh = false;
    cfg.zscore = {};  % ? z-scoring may obscure important differences in VP presence!
%     cfg.zscore_win = [-0.2 0.6];
%     cfg.zscore_cond  = true;
    [data,tbl_depths] = VPmonkey_mergeSesh(cfg);
    length(data.LFP)

end