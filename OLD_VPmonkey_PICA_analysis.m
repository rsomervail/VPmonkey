%
%
%       super rough version using the lw version, should do in EEGlab or something
%       
%
%%
clc 
clearvars -except lwdata

d = sqz(lwdata.data);
h = lwdata.header;

um = h.history.configuration.parameters.ICA_um;
mm = h.history.configuration.parameters.ICA_mm;

nchans = size(um,2);
ncomps = size(um,1);
nconds = 3;



comps_long = [];
for c = 1:nconds

    trial = sqz( d(c,:,:) );

    comps(:,:,c) = um * trial;
    comps_long = [comps_long, um * trial];

end
nsamps = size(comps,2);

contribs_tc = abs(comps) ./ sum(abs(comps), 1); % get normalized contribution timecourse
[maxVals, maxInd] = max(contribs_tc,[],3);
for k = 1 :ncomps

    % compute selectivity metric
    selectivity_tc(k,:) = ( (maxVals(k,:)*3) -1) /2; % ranges from 0 (no selectivity) to 1 (one condition totally dominates)

    % compute mean power of component across time
    temp = mm(:,k) * comps_long(k,:);
    temp = reshape(temp, [nchans, nsamps, nconds]);
    meanposvar_tc(k,:,:) = sum( abs( temp ) ); % ? could also use vector magnitude (norm(x)), might give different results
    
end

% compute proportion of total variance across time
totalposvar_tc = sum( meanposvar_tc ) ;
normposvar_tc  = meanposvar_tc ./ totalposvar_tc;


% compute ratios of contribution to each epoch from activations (? i.e. not the usual method from Liang et al 2010 etc)
times = h.xstart : h.xstep : (nsamps-1)*h.xstep + h.xstart;
t1 = findnearest(times,0);
t2 = findnearest(times,0.5);
contribs = sqz( mean(abs(  comps(:,t1:t2,:)  )  , 2));
contribs = contribs ./ sum(contribs,2); % normalize as proportion total contribution
[maxVals, maxInd] = max(contribs,[],2);
for k = 1 :ncomps

    % compute selectivity metric
%     notmaxInd = setdiff([1:3], maxInd(k));
    selectivity(k) = ( (maxVals(k)*3) -1) /2; % ranges from 0 (no selectivity) to 1 (one condition totally dominates)
    
    % compute mean power of component in post-stimulus period
    temp = mm(:,k) * comps_long(k,:);
    temp = reshape(temp, [nchans, nsamps, nconds]);
    meanposvar(k) = sum( abs( temp(:,t1:t2,:) ) , 'all' ); % ? could also use vector magnitude (norm(x)), might give different results#

end  %  figure; plot(meanposvar)


% compute proportion of total variance 
totalposvar = sum( meanposvar );
normposvar  = meanposvar / totalposvar;


% bin proportion of explained signal
[normposvar_bincounts, normposvar_edges, normposvar_bindy] = histcounts(normposvar);
% compute mean selectivity score by each bin
normposvar_bin_sel = nan(length(normposvar_bincounts),1);
for b = 1:length(normposvar_bincounts)
    if any(normposvar_bindy==b)
        normposvar_bin_sel(b) = mean(  selectivity(normposvar_bindy==b)   );
    end
end

% bin selectivity scores
binWidth = 0.05; % on a scale from 0 to 1
[sel_bincounts, sel_edges, sel_bindy] = histcounts(selectivity, 'BinLimits', [0 1], 'BinWidth',binWidth);
% compute ERP total variance explained by each selectivity bin
sel_bin_normposvar = zeros(length(sel_bincounts),1);
for b = 1:length(sel_bincounts)
    if any(sel_bindy==b)
        sel_bin_normposvar(b) = sum(  normposvar(sel_bindy==b)   );
    end
end

% compute cumulative explained variance by selectivity (like JNP 2009 but showing a whole range of thresholds)
for b = 1:length(sel_bincounts)
    sel_bin_normposvar_cumulative(b) = sum(sel_bin_normposvar(b:end));
end


%% plot stuff 
% plot proportion of explained signal by bin
figure; 
subplot(2,1,1); histogram(normposvar, 'BinEdges',normposvar_edges); ylabel('# components'); xlabel('% explained signal')
binMids = mean(diff(normposvar_edges))/2 : mean(diff(normposvar_edges)) : normposvar_edges(end);
subplot(2,1,2); bar( binMids, normposvar_bin_sel, 1  ); ylim([0, 1]); xlabel('% explained signal'); ylabel('mean selectivity score')

% plot selectivity bins against component count and proportion of signal explained
figure; 
subplot(2,1,1); histogram( selectivity , 'BinWidth', binWidth); xlim([0, 1]); xlabel('selectivity score (0 to 1)'); ylabel('# components')
binMids = binWidth/2 : binWidth : 1;
subplot(2,1,2); bar( binMids, 100*sel_bin_normposvar, 1  ); xlim([0, 1]); xlabel('selectivity score (0 to 1)'); ylabel('total % signal explained')

% plot cumulative explained variance by selectivity  % ! draw area too
figure;
binMids = binWidth/2 : binWidth : 1;
plot(binMids, sel_bin_normposvar_cumulative)

% plot selectivity scores against component variance
figure; scatter(  normposvar , selectivity); 
ax = gca; ax.XAxis.Scale = 'log'; lsline
[r ,p] =  corr( normposvar', selectivity', 'type','Spearman');
title(sprintf('r = %.2f   p = %.6f',r, p))

% % plot selectivity scores against component variance (removing the components explaining small amounts of variance)
% thresh = 10/100; % choose percentage of largest components (?! probs better to just use an absolute "stim-evoked" threshold like in JNP 2009)
% [~, temp] = sort(normposvar,'descend'); temp = temp( 1 : ceil(thresh*numIC) );
% normposvar_thresh  = normposvar(temp);   
% selectivity_thresh = selectivity(temp); 
% % figure; subplot(2,1,1); histogram(normposvar_thresh); subplot(2,1,2); histogram(selectivity_thresh);
% figure; scatter(  normposvar_thresh , selectivity_thresh); 
% ax = gca; ax.XAxis.Scale = 'log'; lsline
% [r ,p] =  corr( normposvar_thresh', selectivity_thresh', 'type','Spearman');
% title(sprintf('r = %.2f   p = %.6f',r, p))

% plot piecharts of normalized contributions for the top N components
topNcomps = ncomps;
figure;
for k = 1:topNcomps
    subplot(4,6,k)
    pie(contribs(k,:))
    title(['c' num2str(k) ' sel = ' num2str(selectivity(k),2)])
end



