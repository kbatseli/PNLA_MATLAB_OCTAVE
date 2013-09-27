function polysys = roots2sys(roots,varargin)
% polysys = roots2sys(roots,multisys)
% -----------------------------------
%
% Finds a multivariate, polynomial system polysys that vanishes at 
% the given roots with given multiplicities in multisys.
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations
%
% roots     =   matrix, each row corresponds with a root, each column
%               corresponds with an indeterminate
%
% multisys  =   cell, polysys-like cell that contains for each root the
%               multiplicity structure, optional.
%
% CALLS
% -----
%
% getAsys.m, reduceAsys.m, fite.m, makeRoot.m, intRound.m
%
% Kim Batselier, 2011-08-04

[d n] = size(roots);

K = makeRoot(d,roots);

LK = null(K')';

R = intRound(flipud((fliplr(rref(fliplr(LK))))));

for i = 1 : size(R,1),A(i,find(R(i,:),1,'last')) = 1;end

Asys =reduceAsys(getAsys(A,n));

[h c indices] = hasAllPureComponents(Asys);

polysys = cell(size(Asys,1),2);
    
for i = 1 : size(Asys,1)
    for j = 1 : size(R,1)
            if sum(frte(find(R(j,:),1,'last'),n) == Asys{i,2}) == n
                % we found the polynomial in R that has the right LT
                coefI = find(R(j,:));
                polysys{i,1} = R(j,coefI);
                for k = 1 : length(coefI)
                    polysys{i,2} = [polysys{i,2} ; fite(coefI(k),n)];
                end
                break
            end
        end
    end
    
%end
