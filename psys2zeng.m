function zeng = psys2zeng(psys)
% zeng = psys2zeng(psys)
% ----------------------
%
% Convert a polysys representation of a multivariate polynomial into a zeng
% matrix.
%
% zeng  =   matrix, each column corresponds with a term, first n rows are
%           exponent, last row is coefficient
%
% psys      =   cell, polysys representation of the polynomial p
%
% CALLS
% -----
%
% Kim Batselier, 2012-02-16

n_terms = length(psys{1,1});
n = size(psys{1,2},2);

zeng = zeros(n+1,n_terms);

for i = 1 : n_terms
    zeng(1:n,i) = psys{1,2}(i,:);
    zeng(n+1,i) = psys{1,1}(i);
end

end