function nz = nnzM(polysys,d)
% nz = nnzM(polysys,d)
% --------------------
%
% This function returns the total number of nonzero elements of M(d),
% constructed from polysys.
%
% nz        =   scalar, total number of nonzero elements of M(d)
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d         =   scalar, desired maximum total degree of matrix M
%
% CALLS
% -----
%
% Kim Batselier

n=size(polysys{1,2},2);
neq=size(polysys,1);
nz=0;
dorig=zeros(1,neq);

for i=1:neq
    dorig(i) = max(sum(polysys{i,2},2));
    nz=nz+length(polysys{i,1})*nchoosek(d-dorig(i)+n,n);
end

end