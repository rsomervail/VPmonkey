

%
%       rickoplot_MK( vals, s) 
% 
%       monkey-optimised version of my usual EEG topoplotting wrapper
% 
%       s.clims 
%       s.colormap
%       s.lay
%       s.gridscale 
%
%
%%  

function  rickoplot_MK( vals, s) 

% takes a vector of data and plots the topo in high res
%%%%%%%%%%%%
% settings %
clims = s.clims;
colormap = s.colormap;
layoutfile = s.lay;

if ~isfield(s,'gridscale')
    s.gridscale = 700;
end
    
if isfield(s,'electrodes')
   electrodes = s.electrodes;
else
   electrodes = 'on';
end


if ischar(layoutfile)
    temp = load(layoutfile);
    if isfield(temp,'lay')
        layout = temp.lay;
    end
    if isfield(temp,'EEG')
        layout = temp.EEG.chanlocs;
    end
    if isfield(temp,'chanlocs')
        layout = temp.chanlocs;
    end
else
    layout = layoutfile; % if chanlocs
end

markerSize  = 20;
lineWidth   =  4; % ears & nose
patchWidth  = 3.7; % scalp radius
contWidth   = 1.5; % contour lines
%%%%%%%%%%%%

 
    rs_topoplot( vals, layout ,...
           'maplimits', clims... 'absmax'
           , 'style', 'both'... % both = color + contour
           ,'shading', 'flat' ...
           , 'colormap',  colormap ...
           ,'numcontour', s.numcontours ...
           ,'gridscale', s.gridscale ... % default = 67, I've been using 700..
           ,'emarker', {'.','k',markerSize,1}  ... % third variable is marker size
           ,'verbose', 'off' ...
           ,'headrad' , 0 ...
           ,'electrodes', electrodes ...
              );
%            ,'pmask',   [ones(1,size(d,1))] ... % set some channels to zero and some to 1 
 

    % change ylims to reveal nose
    xlim( [-0.5373 0.5373] )
    ylim( [-0.5835 0.6631] )
    
    % change line thickness
    if s.numcontours > 0
        a = gca;
        a = a.Children;    
        a(3).LineWidth = contWidth; % contour width
    end
    
    
end 