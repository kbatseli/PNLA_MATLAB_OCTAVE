function [israd e] = isRadical(polysys)
% [israd e] = isRadical(polysys)
% ------------------------------
% Returns true when the given polynomial system is radical. Only makes
% sense if the polynomial Ideal is zero-dimensional.
%
% israd     =   boolean, 1 if polysys is radical, 0 if not
%
% e         =   vector, each row contains 2-norms of LCM and GCD residuals
%               respectively
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
% CALLS
% -----
% 
% getLCM.m, mconv.m, intRound.m, prune.m, diffPoly.m, elim.m
%
% Kim Batselier, 2011-08-02

% number of indeterminates
n = size(polysys{1,2},2);

israd = 1;

e=zeros(n,2);

for i = 1 : n
    
    [punivar angle tol] = univarpol(polysys,i);
    tol=tol(end);
    % then make the monic generators square-free
    % in order to avoid h having mixed terms we convert everything to the
    % univariate case
    punivar=punivar(1:find(abs(punivar)>tol,1,'last')); % prune numerical zero higher order terms
    punivar(abs(punivar)<tol)=0; % set remaining numerical zeros to exact zero
    ppunivar = prune(diffPoly(punivar,1,1),1);
    ppunivar = ppunivar/norm(ppunivar);
    [lcm h e(i,1)] = getLCM(vec2polysys(ppunivar,1),vec2polysys(punivar,1));
    h=h/norm(h);
    
    % computing the GCD-residual norm gives us an indication of how sure we
    % are about the correctness of LCM and h
    gcd = ppunivar/getM(vec2polysys(h,1),length(ppunivar)-1);
    e(i,2) = norm(ppunivar-gcd*getM(vec2polysys(h,1),length(ppunivar)-1));
        
%     % debug function to have some idea on the conditioning of the problem
%     if  (size(lcm,1) == size(prod,1))
%         abs(prod'*lcm)-1
%     end
    
    % either lcm and p*pp have different size or their difference differs
    % from zero within a tolerance of 1e-10    
%     if  (length(lcm)~= length(pmpp)) || abs(pmpp*lcm')-1 > tol
    if  (length(h)~= length(ppunivar)) || abs(ppunivar*h')-1 > tol
        israd = 0;
        return
    end
        
end
