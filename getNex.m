function nex = getNex(polysys,d)
% nex = getNex(polysys,d)
% -----------------------
% Returns the number of indeterminates (exponents) in each syzygy of the
% left null space of the kernel K.
%
% nex       =   vector, each row corresponds with a syzygy of the left null
%               space of the kernel K and holds the number of
%               indeterminates that appear in the syzygy
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d         =   scalar, degree for which Macaulay matrix M is made
%
% CALLS
% -----
% 
% getM.m, fite.m, intRound.m
%
% Kim Batselier, 2011-08-02

n = size(polysys{1,2},2);
M = getM(polysys,d);
K = null(M);
LK = intRound(null(K'))';

R = fliplr(intRound(rref(fliplr(LK))));

expused = zeros(size(R,1),n);

for i = 1 : size(R,1)
    coefI = find(R(i,:));
    for j = 1 : length(coefI)
        expused(i,:) = expused(i,:) + fite(coefI(j),n);
    end
    nex(i,1) = sum(expused(i,:) > 0);
end

end
