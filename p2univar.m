function punivar = p2univar(p,n,var)
% punivar = p2univar(p,n,var)
% ---------------------------
% Converts a univariate polynomial from a vector representation in n
% indeterminates to a vector representation in 1 indeterminate.
%
% p         =   vector, univariate polynomial in basis of n indeterminates
%
% n         =   scalar, number of indeterminates
%
% var       =   scalar, contains index of variable of p
%
% CALLS
% -----
% 
% feti.m, deg.m
%
% Kim Batselier, 2011-07-27
d = deg(p,n);

if issparse(p)
    punivar = spalloc(1,d+1,d+1);
else
    punivar = zeros(1,d+1);
end

for i = 0 : d
    punivar(i+1) = p(feti([zeros(1,var-1) i zeros(1,n-var)]));
end


end
