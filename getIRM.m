function [Mb Q e] = getIRM(polysys,d)
% [Mb e] = getIRM(polysys,d)
% --------------------------
%
% Uses a QR decomposition to determine linear independent rows of Md of a
% given polynomial system polysys and degree d
%
% Mb        =   matrix, contains only the rows of getM(polysys,d) which
%               span the row space, all other spurious rows are deleted
%
% Q         =   matrix, an orthogonal basis for row space of
%               getM(polysys,d)
%
% e         =   vector, contains the indices of the rows of getM(polysys,d)
%               which are used for Mb
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d         =   scalar, desired maximum total degree of matrix M
%
% CALLS
% -----
%
% getM.m
%
% Kim Batselier, 2011-07-10


M = getM(polysys,d);
[Q R E] = qr(M',0);

r = abs(diag(R));
tol = max(size(R)) * eps(max(r));
rankM = sum(r > tol);

e = sort(E(1:rankM));

Mb = M(e,:);

Q = Q(:,1:rankM)';

end

