%
%    fetch EEG peak latencies for all the stimulus modalities for each animal
%       (for use with the scripts including topo-plotting)
%
%       - Richard Somervail, 2025
%%
function topo = VPmonkey_fetchTopoPeaks

    topo.SubM.AUD = [30 95 160]/1000; 
    topo.SubM.SOM = [36 97 177] / 1000;
    topo.SubM.VIS = [22 43  67 98 120 210 290 ]/1000;
    topo.SubT.AUD = [33 79 150]/1000;
    topo.SubT.SOM = [36 120 170] / 1000;
    topo.SubT.VIS = [19 46  110  190 290 ]/1000;

end