


function badtrials_out = VPmonkey_ft_plotTrials(data, ylims, channels, badtrials, layout_file_name)

%% Prepare layout
load(layout_file_name); % this normally load a 'lay' structure
cfg.layout          = lay;
layout              = ft_prepare_layout(cfg, data);

% RS EDIT: REORDER LAYOUT TO MATCH CHANS IN DATA 
newOrder = [ length(layout.label)-2 : -1 : 1     ,  length(layout.label)-1, length(layout.label)  ];
layout.label = layout.label(newOrder);
layout.pos = layout.pos(newOrder,:);

%% convert badtrials to artifact def for visual marking in ft_databrowser
nsamps = size(data.trial{1},2);
indy = find(badtrials);
artdef = [];
for k = 1:length(indy)
    artdef = [artdef;  (indy(k)-1)*nsamps + 1 , (indy(k)-1)*nsamps + nsamps ];
end

%% plot in fieldtrip
    cfg                 = [];
    cfg.layout          = layout;
    cfg.ylim            = ylims;
    cfg.selectmode      = 'markartifact';
    cfg.channel         =  channels;
    cfg.continuous      = 'no';
    cfg.viewmode        = 'vertical';
    cfg.axisfontsize    = 10; 
    cfg.plotlabels      = 'yes';
    cfg.artfctdef.visual.artifact   = artdef;  
    cfg                 = ft_databrowser(cfg,data); 

%% extract bad trials and add to output
for k = 1:length(data.trial)
    trl(k,1) = (k-1)*nsamps + 1;
    trl(k,2) = (k-1)*nsamps + nsamps;
end

% 
artdef =  cfg.artfctdef.visual.artifact;
badtrials_out = false(size(badtrials));
for k = 1 : size(trl,1)  
    if any( artdef(:,1) >= trl(k,1) & artdef(:,2) <= trl(k,2) ) % if trial contains any artifact inside it then reject
        badtrials_out(k) = true;
    end
end





end







