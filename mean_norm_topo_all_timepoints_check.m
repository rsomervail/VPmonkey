
% ? this analysis seems to suggest that EEG timepoints are not excessively weighted in favour of the LFP site at all times
%   ! although note this needs to be run across all sessions in both subs

sizes = size(EEG_all.data);
temp = reshape( EEG_all.data, [sizes(1),sizes(2)*sizes(3)] );
temp = temp ./ vecnorm(temp);
mtopo = mean(temp)';
[~,~,~,stat] = ttest(temp');
ttopo = stat.tstat';

figure;
subplot(1,2,1); rs_topoplot( mtopo, EEG.chanlocs, 'colormap', cmap,   'maplimits', 'absmax'); title('mean')
subplot(1,2,2); rs_topoplot( ttopo, EEG.chanlocs, 'colormap', cmap,   'maplimits', 'absmax'); title('tval')