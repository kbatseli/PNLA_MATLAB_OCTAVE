function Df = diffPoly(f,n,x)
% Df = diffPoly(f,n,x)
% --------------------
%
% Applies the (partial) differential operator on a polynomial f. Graded
% xelicographic ordering is implicitly assumed. Df has the same dimensions
% as f, pruning is not performed.
%
% Df    =   row vector, contains the coefficients of the polynomial
%
%           df
%           --
%           dx
%
% f     =   row vector, contains the coefficients of the multivariate
%           polynomial f
%
% n     =   scalar, number of variables
%
% x     =   scalar, index of the variable to which the differentation needs
%           to take place,
%
%                d        d        d
%           1 = ---, 2 = ---, 3 = ---, etc...
%               dx1      dx2      dx3  
%               
%
% CALLS
% -----
%
% vec2polysys.m
%
% deg.m
%
% getMon.m
%
% Kim Batselier, 2010-01-04

% make sure our f is a column vector
f = f(:);

% check function inputs
if x > n
        error(['You cannot differentiate with respect to x' num2str(x) ', there are only ' num2str(n) ' variables.'])
end

if issparse(f)
    Df = spalloc(length(f),length(f),length(f));
    D = spalloc(length(f),length(f),length(f));
else
    % initialize our differential operator
    Df = zeros(length(f),1);
    D = zeros(length(f),length(f));
end

% determine degree of the polynomial
d = deg(f,n);

% we need a monomial base to construct our operator
base = getMon(d,n);

pos = [];
coeff = [];
for i = 1 : d    
    % for each degree:
    % get the positions of 'x' to the ith degree
    pos = [pos find(base(:,x)==i)'];
    % make a coefficient vector for that degree
    coeff = [coeff i*ones(1,length(find(base(:,x)==i)))];
end

[Y I] = sort(pos);
coeff = coeff(I);

for i = 1 : length(Y)
    D(i,Y(i)) = coeff(i);
end

Df = (D*f)';

end