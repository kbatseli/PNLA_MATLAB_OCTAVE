function [dpi dmi] = getDpi(syzsys,d)
% [dpi dmi] = getDpi(syzsys,d)
% ----------------------------
% 
% Computes for a given syzygy-system the degrees for which positive and
% negative contributions to lcr(d) occur.
%
% dpi       =   vector, contains the degrees for which positive
%               contributions occur
%
% dmi       =   vector, contains the degrees for which negative
%               contributions occur
%
% syzsys    =   cell, polysys object which only contains monomials
%
% d         =   scalar, degree of original polynomial f_i
%
% CALLS
% -----
%
% Kim Batselier, 2012

n = size(syzsys{1,2},2);
n_mon = size(syzsys,1);

expmatrix = zeros(n_mon,n);
% convert the exponent cells of syzsys into a matrix
for i = 1 : n_mon
    expmatrix(i,:) = syzsys{i,2};
end

% initialize the output arguments
dpi = d*ones(1,n_mon+2^(n_mon-1)-n_mon);
dmi = d*ones(1,2^(n_mon-1)-1);

% first we count the basis positive contributions
for i = 1 : n_mon
    dpi(i) = dpi(i)+sum(syzsys{i,2});
end

dpicounter = n_mon+1;
dmicounter = 1;

% initialize combinations
combinations = zeros(2^n_mon-1,n_mon);
degs = zeros(1,2^n_mon-1);

for i = 1:2^n_mon-1    % need all binary strings from 1 to 2^n_mon-1
    binstring = dec2bin(i,n_mon);
    combinations(i,strfind(binstring,'1')) = 1;
    degs(1,i) = sum(combinations(i,:));
end

combinations = logical(combinations); % convert to logical indices

for i = 2:n_mon
    indices = find(degs==i);
    for j = 1:length(indices)
        lcm = max(expmatrix(combinations(indices(j),:),:));
        if ~mod(i,2)
            dmi(dmicounter) = sum(lcm)+d;
            dmicounter = dmicounter + 1; 
        else
            dpi(dpicounter) = sum(lcm)+d;
            dpicounter = dpicounter + 1;            
        end
    end
    
end

dpi = sort(dpi);
dmi = sort(dmi);
end
