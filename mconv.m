function y = mconv(p1,p2,n)
% function y = mconv(p1,p2,n)
% ---------------------------
%
% Generalization of the convolution algorithm to multiply multivariate
% polynomials with each other.
%
% y     =   row vector, contains the coefficients of the multivariate
%           polynomial p1*p2, default: nsgrlex monomial ordering
%
% p1    =   row vector, contains the coefficients of the multivariate
%           polynomial p1, default nsgrlex monomial ordering
%
% p2    =   see p1
%
% n     =   scalar, number of variables
%
% CALLS
% -----
%
% vec2polysys.m
%
% getMDim.m
%
% getM.m
%
% Kim Batselier, 2009-11-17,

% p1 and p2 need to be represented in a full canonical base, so their
% length should be nchoosek(n+d,d)

% make sure p1 and p2 are row vectors
p1 = p1(:)';
p2 = p2(:)';

% convert p1 to polysys
p1sys = vec2polysys(p1,n);

% % number of rows of M needs to be equal to number of columns of p2
% l2 = length(p2);
% d = 0;
% while getMDim(p1sys,d) ~= l2
%     d = d+1;
% end

if issparse(p1) & issparse(p2)
    
    M = getM(p1sys,deg(p1,n)+deg(p2,n),1,p2);
else
    
    M= getM(p1sys,deg(p1,n)+deg(p2,n),0,p2);
end

y = p2*M;

end