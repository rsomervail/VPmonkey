%
%   load intracranial data no filter (for spikes) and import to EEGLAB
%
%%
function [outM, outT] = VPmonkey_loadSPK( path , set )

conds       = set.conds;
trignums    = set.lfp_trig;

%% load data and info
data = TDTbin2mat(path);
d    = data.streams.RAWs.data;
fs   = data.streams.RAWs.fs;
ents.codes = data.epocs.Trig.data;
ents.lats  = data.epocs.Trig.onset;

%% load into EEGlab 
assignin( 'base', 'd', d );
EEG = pop_importdata('data', 'd', 'dataformat', 'array', ...
    'nbchan', 10, 'pnts', size(d,2), 'srate', fs);
evalin('base', 'clear d') % clear d variable from global workspace

% import events to eeglab 
for e = 1:length(ents.codes)
    ents.types{e,1} = conds{trignums==ents.codes(e)};
    [~, ents.latspnts{e}] =  min( abs( ents.lats(e) - (EEG.times/1000) ));
end
EEG = eeg_addnewevents(EEG, ents.latspnts, ents.types );

% loop through events and interpolate over shock artifacts
for e = 1:length(EEG.event)
    if strcmp('SOM', EEG.event(e).type)
        indy = EEG.event(e).latency;
        indys = indy+15:indy+30;
        for c = 1:EEG.nbchan
%             figure; plot(  EEG.times(indy-100:indy+100)/1000 ,  EEG.data(c,  indy-100:indy+100  )  ); hold on
            EEG.data(c,indys) = linspace( EEG.data(c,indys(1)) , EEG.data(c,indys(end)) ,  length(indys)   );
%                     plot(  EEG.times(indy-100:indy+100)/1000 ,  EEG.data(c,  indy-100:indy+100  )  )
        end
    end
end

%% split channels into M and T
outM = pop_select(EEG, 'channel', [6:10]);
outT = pop_select(EEG, 'channel', [1:5 ]); 


end