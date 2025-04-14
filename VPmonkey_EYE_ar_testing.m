%
%
%       testing the ar_eye signed trial rejection in VPmonkey_mergeSesh
% 
%       1)  ar_eye.thresh_sign = [-0.05,5/100];
%           ar_eye.timewin = [-0.5 3];
%               - still misses quite a few dropouts, probs because they are short, try decreasing %
%               - also misses some dropouts happening at very overall amplitudes, 
%                   might need a derivative/freq approach for these
%               - also misses big increases of pupil diameter that seem weirdly fast and stay up for ages
%                 (these probably aren't physiological - try with (higher) positive thresh as well?)
%                   - I think these are actually the cause of the sharp long-lasting broadening of CI in LME_EYE
% 
% 
%       2)  ar_eye.thresh_sign = [-0.05,3/100];
%           ar_eye.timewin = [-0.5 3];
%               - might be removing some visual decreases now ... probably need to exclude that interval
%               - otherwise still doesn't remove some relative dropouts and long very high values but does alright
%               - trying this with LME_EYE_LFP with SubT...
% 
%%
for sb = 1:2

    if sb == 1
        sub = 'SubM';
    else
        sub = 'SubT';
    end

    % ar_EYE
    ar_eye.thresh_sign = [-0.05,3/100];
    ar_eye.timewin = [-0.5 3];
    ar_eye.maxbad  = 0.25;

    cfg = struct; cfg.exportlw = false;
    cfg.sub = sub;
    if exist('ar_eye','var'), cfg.ar_eye = ar_eye; end
    cfg.filetypes = {'EYE'};
    cfg.filt_eye = false; 
    cfg.EYE_include = {'LP_5'}; % only include low-pass filtered pupil data
    cfg.byElectrode = false;
    cfg.average   = false;
    cfg.mergesesh = false;

    [data,tbl_depths] = VPmonkey_mergeSesh(cfg);

end