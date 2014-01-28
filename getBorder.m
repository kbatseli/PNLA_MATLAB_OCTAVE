function border=getBorder(b,n)
% border=getBorder(b,n)
% ---------------------
% Computes the border of a given set of n-variate monomials b.
%
% border    =   vector, contains the indices of the monomials according
%	 			to the degree negative lex (graded xel) monomial ordering
%				that constitue the border of the monomials b,
%
% b         =   vector, contains the indices of monomials according
%	 			to the degree negative lex (graded xel) monomial ordering,
%				
% n         =   scalar, number of variables.
%
% CALLS
% -----
% 
% feti.m

shifts=[];
for j=1:length(b)
    for i=1:n
        shifts = [shifts,feti( fite(b(j),n) + [zeros(1,i-1) 1 zeros(1,n-i)] )];
    end
end
border=setdiff(shifts,intersect(b,shifts));

end
