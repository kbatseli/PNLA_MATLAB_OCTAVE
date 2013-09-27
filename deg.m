function d = deg(p,n)
% d = deg(p,n)
% ------------
% 
% Returns the degree of the multivariate polynomial p, represented by its
% coefficient vector.
%
% d     =   scalar, the degree of p
%
% p     =   vector, contains the coefficients of the multivariate
%               polynomial.
%
% n     =   scalar, number of variables.
%
% CALLS
% -----
%
% fite.m
%
% Kim Batselier, 2009-12-15

d = sum(fite(find(p,1,'last'),n));


end
