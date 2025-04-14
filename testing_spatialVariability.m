
clearvars -except ALLEEG

figure;
for f = 1:length(ALLEEG)

EEG = ALLEEG(f);

EEG.times = EEG.times*1000; % lw export gives times in seconds instead of ms

%    EEG = pop_select(EEG, 'channel', 1:28);

cfg = [];
cfg.func = @rs_spatialVariability;
% cfg.func = @mean; %  cfg.func_params = {  };
cfg.win = [-0.125 0.125];
% cfg.win = [-0.0977 0.0977];
% cfg.win = [ -0.25  0.25   ]; % ! too wide probs
cfg.fs = EEG.srate;
cfg.toi = round( 1: 0.125*EEG.srate  :EEG.pnts);
% cfg.toi = round( 1: 0.25*EEG.srate  :EEG.pnts);
% cfg.toi = linspace(1, EEG.pnts, 8);
cfg.edge_mode = 'exclude';

times = (EEG.times(cfg.toi))/1000;
out = nan(EEG.trials, length(cfg.toi));
% out = nan(EEG.trials, length(cfg.toi), 28); % for functions giving output per channel
for trl = 1:EEG.trials

    cfg2 = [];
%     cfg.func_params = []; % normalisation

    data = squeeze( EEG.data(:,:,trl))';

    % USE ONLY GFP PEAKS --------------------------------------------------------
    % compute GFP peaks
    gfp = std(data');
    [~,locs] = findpeaks(gfp);
    % set all non-GFP peaks timepoints to nan so they are excluded below 
    indy = true(size(data,1),1);
    indy(locs) = false;
    data(indy, :) = nan;
    % ---------------------------------------------------------------------------

    out(trl,:) = rs_winfun(data, cfg); % run sliding window function

%     trl

end


% baseline correction
bl_indy(1) = min( findnearest(times,  -1.5) );
bl_indy(2) = min( findnearest(times,  -0.5) );
out_bl = out - mean( out(:, bl_indy(1):bl_indy(2)), 2, 'OmitNaN' );

% t-test
[~, pvals, ~, stat] = ttest(out_bl);
tvals_all(f,:) = stat.tstat;

out_bl_mean = mean(out_bl);
out_bl_mean_all(f,:) = out_bl_mean;


% plot(  times, out_bl_mean ); hold on
% figure; 
% subplot(2,1,1); plot(  times, out_bl_mean )
% subplot(2,1,2); plot(  times, pvals ); hold on; plot([xlim],[0.05 0.05])

f
end

% plot all subs
figure; 


out_bl_groupmean = mean( out_bl_mean_all );

% using subject means
% t-test
[~, pvals_group, ~, stat] = ttest(out_bl_mean_all);
tvals_group = stat.tstat;
figure; 
subplot(3,1,1)
for f = 1:length(ALLEEG)
    plot( times, out_bl_mean_all(f,:) ); hold on
end; title('using subject means')
subplot(3,1,2); plot(  times, out_bl_groupmean )
subplot(3,1,3); plot(  times, pvals_group ); hold on; plot([xlim],[0.05 0.05])
% subplot(3,1,3); plot(  times, tvals_group ); % hold on; plot([xlim],[0.05 0.05])

% using subject tvalues (i.e. accounting for consistency of effect in each subject)
% t-test
[~, pvals_group_tvals, ~, stat] = ttest(tvals_all);
tvals_tvals = stat.tstat;
figure; 
subplot(3,1,1)
for f = 1:length(ALLEEG)
    plot( times, tvals_all(f,:) ); hold on
end; title('using subject t-values')
subplot(3,1,2); plot(  times, mean(tvals_all) )
subplot(3,1,3); plot(  times, pvals_group_tvals ); hold on; plot([xlim],[0.05 0.05])
% subplot(3,1,3); plot(  times, tvals_tvals ); % hold on; plot([xlim],[0.05 0.05])


