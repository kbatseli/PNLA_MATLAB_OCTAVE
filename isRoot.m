function [rootOK Z] = isRoot(polysys,root,varargin)
% [rootOK Z] = isRoot(polysys,root,tol)
% -------------------------------------
% Checks whether the given roots are roots of the set of polynomial
% equations in polysys. The root argument has to be a matrix where each row
% is a candidate root.
%
% rootOk    =   column vector, zero entry in row i indicates that root(i,:)
%               was not a root of polysys, 1 indicates it was a root.
%
% Z         =   matrix, each column contains the product of M*v for root(i,:)
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% root      =   matrix, each row consists of a candidate root for polysys.
%
% tol       =   scalar, optional. If norm(Z) < tol then the given candidate is
%               considered to be a root, default: 1e-6.
%
% CALLS
% -----
%
% getM.m, getMon.m
% 
% Kim Batselier, 2009-10-13

if ~isempty(varargin)
    tol = varargin{1};
else
    % when Mv = 0, then right-hand side is a column vector of length size(M,2)
    tol = 1e-6;    % if norm(Mv) bigger than this tolerance then we consider v not a solution of M
end
rootOK = zeros(size(root,1),1);
Z = zeros(size(polysys{1,2},2),size(root,1));

warning off all
[M d] = getM(polysys,0);

monBase = getMon(d,size(polysys{1,2},2));

for j = 1 : size(root,1)                        % for each candidate root
    
    v = zeros(size(monBase,1),1);               % initialize root vector
    for i = 1 : size(monBase,1)                 % iterate over monomial basis
        v(i,:) = prod(root(j,:).^monBase(i,:));
    end
    Z(:,j) = M*v;                               % get the right-hand side
    rootOK(j,1) = norm(Z(:,j)) < tol;       % root-condition
end

warning on all
end
