function [NS normalset] = getNS(polysys,d)
% [NS normalset] = getNS(polysys,d)
% -----------------------------------
%
% Calculates the standard monomials (Normal Set) which are also a basis for
% the quotient space R/I where R refers to the Polynomial Ring and I to an
% Ideal.
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
% getM.m,  
%
% Kim Batselier, 2010-11-04, update: now uses sparseqr and permutations
% as in sparf.m and polyDiv.m

n = size(polysys{1,2},2);

M = getM(polysys,d,1);   

% %% determine rank of Divisor matrix
[Q R P] = qr(M',struct('Q','matrix','permutation','vector'));
rankM = length(find(diag(R)));
V = Q(:,rankM+1:end);

clear Q R P

% Determine normal set that will span the remainder
[Qv Rv Pv] = qr(V',struct('Q','discard','permutation','vector'));
normalset = sort(Pv(1:size(V,2)));
NS = zeros(size(V,2),n);
for i = 1 : size(V,2)
    NS(i,:) = fite(normalset(i),n);
end

end
