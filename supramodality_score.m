



%
function sup = supramodality_score(a)

    %% Option 1 - product of all ratios
%     n = length(a);
%     m = nan(n);
%     for r = 1:n
%         for c = 1:n
%             if r~=c
%                 m(r,c) = a(r) / a(c);
%             end
%         end
%     end
% 
%     up = triu(m); up = up(:)'; up = up(~isnan(up) & up~=0);
%     dn = tril(m); dn = dn(:)'; dn = dn(~isnan(dn) & dn~=0);
%     rats = max([up;dn]);
% 
%     sup = prod(rats);

    %% Option 2 - biggest number divided by mean of the others
%     n = length(a);
%     
%     [m,inds_max] = max(a);
%     inds_other = setdiff(1:n,inds_max);
% 
%     sup = m / mean(a(inds_other));
    
    %% Option 3 - biggest number minus mean of the others
    % (this worked well for my Somervail, 2022 paper)
    n = length(a);
    
    [m,inds_max] = max(a);
    inds_other = setdiff(1:n,inds_max);

    sup = m - mean(a(inds_other));

end




