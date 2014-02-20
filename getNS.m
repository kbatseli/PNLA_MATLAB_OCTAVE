function [NS normalset] = getNS(polysys,d)
% [NS normalset] = getNS(polysys,d)
% -----------------------------------
%
% Obsolete: It is better (more robust) to compute the normal set monomials via
% the canonical decomposition method candecomp.m rather than with this method.
% Calculates the standard monomials (Normal Set) for a given degree d 
% of the Macaulay matrix. These monomials are also a basis for
% the quotient space R/I where R refers to the Polynomial Ring and I to a
% polynomial ideal for a large enough degree d.
%
% NS        =   matrix, each row contains the exponents of a standard
%               monomial
%
% normalset =   vector, contains indices of dependent columns of M
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations that span ideal I.
%
% d         =   scalar, degree for which the normal set needs to be
%               determined
%
% CALLS
% -----
%
% getM.m,  spqr.m (SuiteSparse package)
%
% Kim Batselier

n = size(polysys{1,2},2);

M = getM(polysys,d,1);   

% %% determine rank of Macaulay matrix
[Q R P] = spqr(M',struct('Q','matrix','permutation','vector'));
rankM = length(find(diag(R)));
V = Q(:,rankM+1:end);

clear Q R P

% Determine normal set that will span the remainder
[Qv Rv Pv] = spqr(V',struct('Q','discard','permutation','vector'));
normalset = sort(Pv(1:size(V,2)));
NS = zeros(size(V,2),n);
for i = 1 : size(V,2)
    NS(i,:) = fite(normalset(i),n);
end

end
