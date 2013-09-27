function Asys = getAsys(A,n)
% Asys = getAsys(A,n) or getAsys(a,n)
% -----------------------------------
%
% Returns the polysys cell for the given A matrix containing a monomial
% basis for LT(Id).
%
% Asys      =   cell, polysys cell for the corresponding A matrix
%
% A         =   matrix, monomial basis for LT(Id)
%
% a         =   vector, contains indices of linear independent leading
%               monomials
%
% n         =   scalar, number of unknowns
%
% CALLS
% -----
%
% fite.m
%
% Kim Batselier, 2011-06-20
if sum(A(1,:)) == 1
    Asys = cell(size(A,1),2);
    
    for i = 1 : size(A,1)
        Asys{i,1} = 1;
        Asys{i,2} = fite(find(A(i,:),1),n);
    end
else
    Asys = cell(length(A),2);
   for i = 1:length(A)
       Asys{i,1}=1;
       Asys{i,2}=fite(A(i),n);
   end
end

end
