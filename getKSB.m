function [D I] = getKSB(d,depth,root)
% [D I] = getKSB(d,depth,root)
% ----------------------------
%
% Calculates and returns a Kernel Subspace Basis, consisting of all partial
% derivatives from degree 0 up to degree depth. The maximum degree of the
% monomial base is d. All these vectors of differentiated monomials are
% evaluated in root.
%
% D     =   matrix, each column contains a partial derivative of a monomial
%           basis, evaluated in root
%
% I     =   matrix, each column contains the indices of differentiation of
%           the corresponding column of D, eg. I(:,5) = [1;1] means that
%           D(:,4) = Dxy.
%
% d     =   scalar, degree of the monomial base
%
% depth =   scalar, maximum degree of differentiation
%
% n     =   scalar, number of variables
%
% root  =   vector, point in n-space in which the vectors of D need to be
%           evaluated.
%
% CALLS
% -----
%
% getMon.m
%
% diffBase.m
%
% makeRoot.m
%
% Kim Batselier, 2010-02-17
n = size(root,2);
D(:,1) = makeRoot(d,root);

I = getMon(depth,n);

for i = 2 : size(I,1)
    temp = ones(1,I(i,1));
    for j = 2 : n
        temp = [temp j*ones(1,I(i,j))];
    end        
    D(:,i) = diffBase(d,root,temp);
end

I = I';

end

