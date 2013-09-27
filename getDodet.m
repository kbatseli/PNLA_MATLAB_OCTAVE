function [d p q] = getDodet(polysys)
% [d p q] = getDodet(polysys)
% ---------------------------
% Determines for which degree d the M matrix will be overdetermined. 
%
% d         =   scalar, degree for which M is overdetermined
%
% p         =   scalar, number of rows of overdetermined M
%
% q         =   scalar, number of columns of overdetermined M
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% CALLS
% -----
% 
% Kim Batselier, 2010-10-17

% number of unknowns
n = size(polysys{1,2},2);

% number of equations
n_eq = size(polysys,1);

di = zeros(1,n_eq);
% figure out the degree of each polynomial
for i = 1 : n_eq
    di(i) = max(sum(polysys{i,2},2));
end

d0 = max(di);
d = d0;

% initialize number of rows and columns
p = 0;
for i = 1 : n_eq
    p = p + nchoosek(d-di(i)+n,n);
end
q = nchoosek(d+n,n);

while p <= q
    % update d, p and q
    d = d+1;
    p = 0;
    for i = 1 : n_eq
        p = p + nchoosek(d-di(i)+n,n);
    end
    q = nchoosek(d+n,n);
end

end