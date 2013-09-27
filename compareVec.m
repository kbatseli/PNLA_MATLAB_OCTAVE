function [abse rele] = compareVec(phat,p)
% [abse rele] = compareVec(phat,p)
% --------------------------------
% 
% Compares an estimated polynomial with its exact version.
%
% abse        = scalar, absolute error
%
% rele	      = scalar, relative error
%
% phat        = coefficient vector of multivariate polynomial
%
% p 	      = coefficient vector of multivariate polynomial
%
% REFS
% ----
%
% A geometrical approach to finding multivariate approximate LCMs and GCDs
%
% Kim Batselier, 2012-07-15

pos = find(abs(p)==max(abs(p)));

% scale phat appropiately by dividing by its largest absolute value
phat = p(pos)*phat/phat(pos);

abse = norm(phat-p);
rele = abse/norm(p);

end
