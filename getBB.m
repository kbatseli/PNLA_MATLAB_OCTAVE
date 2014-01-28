function  bb = getBB(polysys,d,b)
% bb = getBB(polysys,d,b)
% -----------------------
% Computes a border pre-basis at degree d for a given polynomial system polysys.
%
% bb	    =   matrix, each row corresponds with a coefficient vector of
%				a polynomial of the border prebasis, 
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations,
%
% d         =   scalar, degree at which border pre-basis needs to be computed,
%
% b         =   vector, contains the indices of monomials according
%	 			to the degree negative lex (graded xel) monomial ordering,			
% CALLS
% -----
% 
% getM.m, getMex.m, getBorder.m, updateN.m

n=size(polysys{1,2},2);
d0=getD0(polysys);
M=getM(polysys,d0);
N=null(M);
% update basis nullspace M(d) recursively
for i=d0+1:d
	N = updateN(N,getMex(polysys,i,i-1));
end

border = getBorder(b,n);
bb = zeros(length(border),length(b)+1);

for i=1:length(border)
	temp = null(N([b border(i)],:)');
	
	if ~isempty(temp)
		if size(temp,2) > 1
			disp(['warning: multiple solutions detected for leading term: ' num2str(border(i))])
			bb(i,:) = temp(:,1);
		else
			bb(i,:) = temp;
		end
	end

end
