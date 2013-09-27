function y = prune(p,n)
% y = prune(p,n)
% --------------
% 
% Prunes zeros away from a multivariate polynomial coefficient vector such
% that the vector has a minimum length.
%
% y         =   row vector, vector p with unnecessary zeros removed.
%
% p         =   vector, multivariate polynomial coefficient vector
%
% n         =   scalar, number of variables%
%
% CALLS
% -----
%
% getMon.m
%
% Kim Batselier, 2009-12-09; updated 2011-08-03: returns a row vector now

y = p(1:nchoosek( sum(fite(find(p,1,'last'),n))+n,n  ));

end
