function [dp dm] = mergeDpi(dpcell,dmcell)
% [dp dm] = mergeDpi(dpcell,dmcell)
% ---------------------------------
% 
% Merges the dpi and dmi's of different syzygy-systems in order to obtain
% the full dp and dm as found by aln.m
%
% dp        =   vector, contains the degrees for which positive
%               contributions occur
%
% dm        =   vector, contains the degrees for which negative
%               contributions occur
%
% dpcell    =   cell, each entry contains the dpi of i'th syzygy-system
%
% dmcell    =   cell, each entry contains the dmi of i'th syzygy-system
%
% CALLS
% -----
%
% Kim Batselier, 2012

dpi =[];
dmi =[];

dp = [];
dm=[];

for i = 1 : size(dpcell,2)
   dpi = [dpi,dpcell{i}];
   dmi = [dmi,dmcell{i}];   
end

if isempty(dmi)
    dmin = min(dpi);
    dmax = max(dpi);
else
    dmin = min(min(dpi),min(dmi));
    dmax = max(max(dpi),max(dmi));
end

for i = dmin:dmax
   if sum(dpi==i)-sum(dmi==i) < 0
       dm = [dm i*ones(1,sum(dmi==i)-sum(dpi==i))];
   elseif sum(dpi==i)-sum(dmi==i) > 0
       dp = [dp i*ones(1,sum(dpi==i)-sum(dmi==i))];
   end
end




end