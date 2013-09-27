function [dp dm] = aln(polysys,dmax)
% [dp dm] = aln(polysys,dmax)
% ---------------------------
% Analyze Left Nullspace of Macaulay matrix. Calculates the Left Hilbert
% Polynomial for a given polynomial system.
%
% dp        =   vector, contains degrees that contribute with + terms to
%               the Left Hilbert Polynomial
%
% dm        =   vector, contains degrees that contribute with - terms to
%               the Left Hilbert Polynomial
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
% getD0.m, getMDim.m, getM.m, getMex.m
% 
% Kim Batselier, 2011-06-09, update 2013-03-25: uses updateN now to
% determine r(d)

d = getD0(polysys);
dp = [];
dm = [];
n = size(polysys{1,2},2);

[p q] = getMDim(polysys,dmax);

[M ] = getM(polysys,0);

for i = 1 : dmax % new algorithm - we always need to start from degree 1
    [p q] = getMDim(polysys,i);
    temp = M;
    [Mex ] = getMex(polysys,i,i-1);
    M = zeros( p,q );
    M(1:size(temp,1),1:size(temp,2)) = temp;
    M(size(temp,1)+1:end,:) = Mex;
%     Ms = [Ms zeros(size(Ms,1),size(Mex,2)-size(Ms,2));Mex];    
    clear Mex
%     M = getM(polysys,i);
    if ~isempty(M)
        lcr = size(M,1)-rank(M);
        lcrhat = evalLcr();
        e = lcrhat-lcr;
        if e > 0       
            for k = 1 : abs(e)
                dm = [dm i];
            end
        elseif e < 0
            for k = 1 : abs(e)
                dp = [dp i];
            end
        end
    end
%     [length(dp) length(dm) length(dp)-length(dm)]
end

    function lcrhat = evalLcr
        lcrhat = 0;
        for j = 1 : length(dp)
            lcrhat = lcrhat + nchoosek(i-dp(j)+n,n);
        end
        for j = 1 : length(dm)
            lcrhat = lcrhat - nchoosek(i-dm(j)+n,n);
        end        
    end

end
