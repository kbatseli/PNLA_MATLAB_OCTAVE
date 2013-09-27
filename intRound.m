function y = intRound(x,varargin)
% y = intRound(x) or y = intRound(x,tol)
% --------------------------------------
% this function rounds the numbers of in x to the nearest integer when the
% difference between the number and the nearest integer is less than a
% certain tolerance (default 1e-10).
%
% y     =   matrix, copy of x, rounded off to nearest integer if applicable
%
% x     =   matrix
%
% tol   =   scalar, if difference is smaller than this number then the
%           difference is considered to be zero, default: 1e-10
%
% CALLS
% -----
%
% Kim Batselier, 2010-01

if isempty(varargin)
    zerotol = 1e-10;
else
    zerotol = varargin{1};
end

y=(abs(x-round(x)) <= zerotol).*round(x) + (abs(x-round(x)) > zerotol).*x;


end