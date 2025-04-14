%
%
%
%   VPmonkey check mean modality offset
%       - takes lw files with only the relevant modality and one channel reflecting the intended metric, already z-scored appropriately
%
%
%
%%

d = sqz(lwdata.data);
h = lwdata.header;

%!!! didn't bother doing this yet, just used rough estimate