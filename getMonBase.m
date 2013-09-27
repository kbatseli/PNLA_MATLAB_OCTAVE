function base = getMonBase(d,n)
% base = getMonBase(d,n)
% ---------------------
%
% Returns a set of monomials of total degree d in n variables. Each row of
% mon refers to a n-tuple of exponents of a monomial whereby each column
% corresponds to a variable.
%
%
% base      =   matrix of lexicographic ordered monomials of degree d and
%               n exponents
%
% d         =   scalar, maximum total degree of the monomials in the base
%
% n         =   scalar, number of exponents
%
% example: d= 2, n = 3
%
% base =
% 
%      2     0     0
%      1     1     0
%      1     0     1
%      0     2     0
%      0     1     1
%      0     0     2
%
% CALLS
% -----
%
% getMonBase.m
%
% Kim Batselier 2009-10

% if n == 1
%     base = d;
% else
%     base = [];
%     for i = d :-1: 0
% %         temp = getMonBase(d-i,n-1);
% %         base = [base; i*ones(size(temp,1),1) getMonBase(d-i,n-1)];
%         base = [base; i*ones(nchoosek(d-i+n-2,n-2),1) getMonBase(d-i,n-1)];
%     end
% end
% disp(['d is ' num2str(d) ', n is ' num2str(n)])
if n == 1
    base = d;
else
    base = [d zeros(1,n-1)]    ;
    for i = d-1 :-1: 0
%         temp = getMonBase(d-i,n-1);
%         base = [base; i*ones(size(temp,1),1) getMonBase(d-i,n-1)];
        base = [base; i*ones(nchoosek(d-i+n-2,n-2),1) getMonBase(d-i,n-1)];
    end
end

end