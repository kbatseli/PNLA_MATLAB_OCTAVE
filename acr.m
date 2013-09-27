function [crp crm] = acr(polysys,dmax)
% [crp crm] = acr(polysys,dmax)
% -----------------------------
% Analyze Co-Rank of Macaulay matrix. Calculates the Right Hilbert
% Polynomial for a given polynomial system.
%
% crp        =  vector, contains degrees that contribute with + terms to
%               the Right Hilbert Polynomial
%
% crm        =   vector, contains degrees that contribute with - terms to
%               the Right Hilbert Polynomial
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% dmax      =   scalar, desired maximum total degree of matrix M
%
%
% CALLS
% -----
%
% getD0.m, getMDim.m, getM.m
% 
% Kim Batselier, 2011-06-09

d = getD0(polysys);
crp = [];
crm = [];
n = size(polysys{1,2},2);

[p q] = getMDim(polysys,dmax);

for i = 1 : dmax % new algorithm - we always need to start from degree 1
    
    M = getM(polysys,i);
    if ~isempty(M)
        cr = size(M,2)-rank(M);
        crhat = evalCr();
        e = crhat-cr;
        if e > 0       
            for k = 1 : abs(e)
                crm = [crm i];
            end
        elseif e < 0
            for k = 1 : abs(e)
                crp = [crp i];
            end
        end
    end
    
end

    function crhat = evalCr
        crhat = 0;
        for j = 1 : length(crp)
            crhat = crhat + nchoosek(i-crp(j)+n,n);
        end
        for j = 1 : length(crm)
            crhat = crhat - nchoosek(i-crm(j)+n,n);
        end        
    end

end
