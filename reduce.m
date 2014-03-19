function ar = reduce(a,n)
% ar = reduce(a,n)
% ----------------
%
% Reduces a given set of monomials to its smallest subset such that each
% monomial in a is divisible by at least one element of ar.
%
% ar    =   vector of indices corresponding with monomials degree negative
%           lex (graded xel) ordered
%
% a     =   vector of indices corresponding with monomials degree negative
%           lex (graded xel) ordered
%
% n     =   scalar, number of variables
%
% CALLS
% -----
%
% getExpM.m, feti.m
%
% Kim Batselier, 2014-02-28

Am=getExpMat(a,n);
ar=[];

counter = 1;
while ~isempty(Am)

    ar(counter)=feti(Am(1,:));
    
    verschil = Am-ones(size(Am,1),1)*Am(1,:);
        
    % remove derivatives of baseA{1} from Am
    indices = verschil(:,1) >= 0;
    for i = 2 : n
        indices = indices & (verschil(:,i)>=0);
    end
    Am(indices,:) = [];
    counter = counter + 1;

end

end