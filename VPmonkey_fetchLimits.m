% 
% 
%       global function for fetching a structure with appropriate x and y limits for each datatype
% 
%           assumes no z-scoring on the estimated means to be plotted
%           but z-scoring of predictors should be fine
% 
%       - Richard Somervail, 2025
%%
function lim = VPmonkey_fetchLimits

% xlims
lim.xlims.plot.EEG  = [-0.2 0.6];  
lim.xlims.plot.LFP  = [-0.2 0.6];  
lim.xlims.plot.MUA  = [-0.2 0.6];  
lim.xlims.plot.EYE  = [-0.5   3]; 
lim.xlims.plot.MISC = [-0.2 0.5];
% ylims
% ylims - amplitudes
lim.ylims.plot.EEG  = [-20 20];  
lim.ylims.plot.LFP  = [-140 140];   
lim.ylims.plot.MUA  = [-0.14 0.14];   % ? check later
lim.ylims.plot.EYE  = [-1 1]; 
% topoplot clims - amplitudes  ? considering using rounded absmax now instead
lim.clims.SubM.AUD = [-5 5];
lim.clims.SubM.SOM = [-5 5];
lim.clims.SubM.VIS = [-5 5];
lim.clims.SubT.AUD = [-5 5];
lim.clims.SubT.SOM = [-5 5];
lim.clims.SubT.VIS = [-5 5];



end

