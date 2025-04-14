function SDF = f_calcSDF(spikeTimes, kernelSD, fixedDuration)
% spike times and fixedDuration   in seconds
% kernelSD   in samples (=milliseconds)
%
% ouptut = SDF, sampled at 1000Hz
%
% v01 20140513
% v02 20140715


if nargin==1,
    kernelSD = 30;  % default kernel width: 30 ms
    fixedDuration = 0;
elseif nargin==2,
    fixedDuration = 0;
end

    % init sdf
    npts = kernelSD*4;
    basicGauss = normpdf(-npts:npts, 0, kernelSD);
    
	% generate output
	spikesT = round(spikeTimes*1000);   % spike times, in ms
	
	if fixedDuration>0,
		S = round(fixedDuration*1000)+1;
	else
		S = max(spikesT)+npts;
	end
	
	sdfArr = zeros(1,S);
	if ~isempty(spikesT),
        
        sdfArr(1,floor(spikesT)+1) = 1;
        sdf_ = 1000*conv(sdfArr, basicGauss)/sum(basicGauss);
        SDF = sdf_(npts+1:1:end-npts);
        
    else
        
		SDF = zeros(1,S);
        
	end
	
	
end