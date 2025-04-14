function [ data_out, interpolated_channels ] = RS_Giac_EEGinterpolation( data, channels, range, layout_file_name, neigh_distance, plotData )
%%
% Function for interpolation of EEG electrodes
% 'channels' specifies which channels have to be selected for the
% interpolation. 'neigh_distance' specifies how close electrodes have to be
% in order to be included in the interpolation. 
% The functions works recursively. You can gradually include more channels
% to be interpolated. Note that, on every instance, the function will
% implement the interpolation ex-novo (including all channels selected so far) to avoid leaking signal from
% channels-to-be-interpolated into clean ones.
% 
% 
%       RS edit: various edits: 
%           - added option to skip data plotting steps
%           - had issue where lay file and locs file are reversed order compared to my data, so this checks that match with the actual conve
%           - updated GUI to make more user friendly
% 
%
%% Giacomo Novembre
 
chan_str    = data.label;

%% Prepare layout
load(layout_file_name); % this normally load a 'lay' structure
cfg.layout          = lay;
layout              = ft_prepare_layout(cfg, data);

% RS EDIT: REORDER LAYOUT TO MATCH CHANS IN DATA 
newOrder = [ length(layout.label)-2 : -1 : 1     ,  length(layout.label)-1, length(layout.label)  ];
layout.label = layout.label(newOrder);
layout.pos = layout.pos(newOrder,:);


%% Provide tips about suspicious electrodes... based on Automatic identification of NOISY electrodes
[ out_ch_noisy ]    = Giac_EEG_CatchNoisyElectrodes( data, channels, 3, 'recursive' ); % find noisy electrodes
badchanindy = cellfun( @(x) find(strcmp(data.label, x)) , out_ch_noisy)'; % get indices of these channels for later plotting and conditionals
fprintf( ['GIAC: suspicious electrodes: ' strjoin(out_ch_noisy, ' ') '\n'] ); % RS edit: fixing some bug with conversion from cell to char (likely due to older matlab version)

%% Loop in case of multiple checks... 
loop                  = 0;
interpolated_channels = [];
original_data         = data; % this should never be changed, it is the data as inputted (no change is implemented)
tmp_data              = data; % this is the data you can change to evaluate the effect of the interpolation

while loop < 1

    %% Start with plotting
    if plotData
        cfg                 = [];
        cfg.layout          = layout;
        cfg.ylim            = range;
        cfg.continuous      = 'no';
        cfg.selectmode      = 'markartifact';
        cfg.channel         = channels;
        cfg.viewmode        = 'vertical';
        cfg.axisfontsize    = 10; 
        cfg.plotlabels      = 'yes';
        cfg.artifact        = [];  
        % RS EDIT: made it so that bad channels are marked in red, others are blue
        cfg.linecolor = zeros(length(data.label),3);
        for c = 1:length(data.label)
            if ismember(c, badchanindy)
                cfg.linecolor(c,1) =  1; % red
            else
                cfg.linecolor(c,3) =  1; % blue
            end
        end
        %%%%%%
        cfg = ft_databrowser(cfg,tmp_data); % CFG contains artifacts (previous one is for plotting only)
    end

    %% Question about interpolation
    if ~isempty(interpolated_channels)
        display(['GIAC: currently you are interpolating electrode : ' strjoin(interpolated_channels', ' ') ]);
    else
        disp 'GIAC: currently interpolating no electrodes'
    end
    if plotData
        Question            = listdlg('PromptString','Channels to interpolate?','SelectionMode','multiple','ListString',{'yes','no'});
    elseif ~isempty(badchanindy)
        Question = 1;
    else
        Question = 2;
    end

    %% Actual interpolation
    if Question==1 && plotData
        ChanInterpol           = listdlg('PromptString','Choose bad channels:','SelectionMode','multiple','ListString',chan_str, ...
            'InitialValue', [badchanindy]  ); % RS EDIT: made so that bad channels auto marked
    elseif Question==1 && ~plotData
        ChanInterpol = badchanindy; % if not plotting data then just assume all reccomended electrodes will be interpolated
    elseif Question==2
        ChanInterpol = [];
    end
    if ~isempty(ChanInterpol)
        interpolated_channels  = [interpolated_channels; chan_str(ChanInterpol)];
                cfg                 = [];
                cfg.layout          = layout;
                cfg.method          = 'distance'; % for prepare_neigh
                cfg.neighbourdist   = neigh_distance;         % results in avg 5 channels
                cfg.neighbours      = ft_prepare_neighbours(cfg, original_data);
                cfg.badchannel      = interpolated_channels; %data.label(ChanInterpol);
                cfg.method          = 'nearest';     
                tmp_data            = ft_channelrepair(cfg, original_data); 
    end

    %% Stop the loop or continue
    if plotData
        Question = listdlg('PromptString','Look again at the data?','SelectionMode','multiple','ListString',{'yes','no'});
    else
        Question = 2;
    end
    if Question == 1 % If another look at the data
        loop = 0;
    elseif Question == 2 % If not... 
        loop = 1;
    end

end % end of loop

data_out = tmp_data;

end

