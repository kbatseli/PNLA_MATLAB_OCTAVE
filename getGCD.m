function [gcd e] = getGCD(psys,qsys,tol)
% [gcd e] = getGCD(psys,qsys,tol)
% -------------------------------
% Returns the approximate greatest common divisor of two given polynomials psys and qsys. The
% mean squared errors (MSE) of the least common multiple and gcd are also
% returned.
%
% gcd       =   vector, coefficient vector of polynomial which is the
%               greatest common divisor of p and q.
%
% e         =   vector, e(1): MSE of LCM, e(2): MSE of GCD
%
% psys      =   cell, polysys representation of the polynomial p
%
% qsys      =   cell, polysys representation of the polynomial q
%
% tol       =   scalar, tolerance: measure of how many digits are not
%               corrputed by noise. Default: eps
%
% CALLS
% -----
% 
% getM.m, getLCM.m, vec2polysys.m, deg.m, getDo.m
%
% Kim Batselier, 2012-01-31

if nargin == 2
    tol = eps;
end

n = size(psys{1,2},2);
[lcm h e(1,1)] = getLCM(psys,qsys,tol);

p = getM(psys,getD0(psys),1)/norm(psys{1,1});

M = getM(vec2polysys(h,n),getD0(psys),1)';
gcd(1,:) = M\p';
e(1,2) = norm(p-gcd(1,:)*M');

end
