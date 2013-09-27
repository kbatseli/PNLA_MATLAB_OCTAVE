function [gsys a b stop dec] = getGsys(polysys,d,varargin)
% [gsys a b stop] = getGsys(polysys,d,tol)
% ----------------------------------------
% Computes a reduced polynomial system such that all leading terms of the
% reduced row echelon form of getM(polysys,d) are divisible by the leading
% terms of gsys. Outdated by reduced canonical decomposition.
%
% gsys      =   cell, polysys cell which contains the reduced polynomial
%               system
%
% a         =   vector, contains indices of reduced linear independent leading
%               monomials
%
% b         =   vector, contains indices of affine linear depedent monomials
%
% stop      =   boolean, is 1 if gsys satisfies the condition of being a
%               Grobner basis of polysys
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d         =   scalar, degree for which the 
%
% tol       =   scalar, optional tolerance for checking numerical zeros.
%               Default: max(size(M))*eps(norm(M)).
%
% CALLS
% -----
%
% fite.m, getM.m, getMon.m, getExpMat.m, hasAllPureComponents.m
%
% Kim Batselier, update 2012-09-18: uses the svd-based reduction now
% instead of the numerical unstable rref.

% first we normalize each polynomial by dividing by its 2-norm, this helps
% for the Cassou system
for i=1:size(polysys,1)
    polysys{i,1}=polysys{i,1}/norm(polysys{i,1});
end

n=size(polysys{1,2},2);
M=getM(polysys,d);
[U S V] = svd(M,'econ');
sig=diag(S);

if nargin == 2
	tol = max(size(M))*eps(sig(1));
else
	tol = varargin{1};
end

checki = [1:size(M,2)]; %contains indices of monomials that need to be checked for linear independence

a=[]; % monomial multiples of elements of a do not need to be checked anymore
b=[];
gsys=[];

% need to take tolerance into account for finding orthogonal basis of the
% row space since the dimension of the row space is the numerical rank
% under the given tolerance
r = sum(sig> max(size(M))*eps(sig(1))); % we assume no higher order terms are introduced by the noise
sig(sig==0)=[];  % remove exact zeros, to avoid division by zero, since anyway they are not counted in the rank-test
[Y I] = max(sig(1:end-1)./sig(2:end)); % compute approxi-rank gaps
if r ~= I  && r ~= size(M,1)  
    disp('Warning: rank estimation is probably wrong')
    disp('------------------------------------------')
    disp(['Estimated rank from approxi-rank gap: ' num2str(I) ' with gap of ' num2str(Y)])
    disp(['Estimated rank from default tolerance: ' num2str(r) ' with gap of ' num2str(sig(r)/sig(r+1))])
    % override tolerance and rank
%     r = I;
%     % tolerance is set to ceil of first nonzero digit of the first 'small'
%     % singular value
%     tol= ceil(sig(r+1)*10^( ceil(-log10(sig(r+1)))))*10^(- ceil(-log10(sig(r+1))));
end
% r = sum(sig>tol);
U= V(:,1:r);
% U = orth(M');
% first iteration is a special case
test=-U*U(1,:)';
test(1,1)=1+test(1,1);
s=norm(test); % returns largest singular value (there is only 1 singular value for first iteration)
if asin(s) < tol
    % '1' lies in the ideal and hence the variety is empty, no need to
    % check the other monomials
    gsys{1,1}=1;
    gsys{1,2}=zeros(1,n);
    return;
else
    b=1;
end

counter =2;
gsyscounter = 1;
while counter <= length(checki)
    indices =[b checki(counter)];
% 	indices=[1:checki(counter)];
% 	% if an A monomial is found, it should be excluded in the indices for the next iteration[
% 	indices(a)=[];

	test=-U*U(indices,:)';
	for j=1:length(indices)
		test(indices(j),j) = 1+test(indices(j),j);
	end
	% compute the sines
	[Y S Z] = svd(test);
	s=diag(S);
	
	%[fite(n,checki(counter)) asin(s(end))]
    dec(counter)=asin(s(end));
	if asin(s(end)) < tol
		a = [a checki(counter)];
		% remove all monomial multiples from checki
		di=sum(fite(checki(counter),n));
		if di<d
			multmon = getMon(d-di,n,1);
			multiples = multmon+ones(size(multmon,1),1)*fite(checki(counter),n);
            for j=1:size(multiples,1)
                checki(checki==fetr(multiples(j,:)))=[];
            end
        end
%         gap = s(end-1)/s(end);        
        p(gsyscounter,indices)=Z(:,end)'; % our required polynomial lies in the kernel
	%[gsyscounter s(end-1) asin(s(end-1)) asin(s(end))]
%         p(gsyscounter,abs(p(gsyscounter,:))<s(end-1)/2)=0; % set coefficients which are smaller than s(end-1)/2 to exact zero!
		gsyscounter = gsyscounter+1;
    else
        b=[b checki(counter)];
    end
	clear test Y S Z s
	counter = counter + 1;

end

% we remove numerical zeros from p according to tol
% p(abs(p)<tol)=0;
gsys = vec2polysys(p,n);

% we need to check whether for all LT's of gsys that the degree of the lcm
% of two elements of gsys does not exceed d
A = getExpMat(a,size(polysys{1,2},2)); % get all exponents of leading monomials
temp= sum(A,2);
if length(temp)==1
    dupper = temp;
else
    dupper = temp(end-1)+temp(end);
end
dmax=0;
setbreak=0;
for i=size(A,1):-1:2
    for j=i-1:-1:1
        if sum(max([A(i,:);A(j,:)])) > dmax
            dmax = sum(max([A(i,:);A(j,:)]));
        end
        if dmax == dupper
           setbreak=1;
           break
        end
    end
    if setbreak
        break
    end
end
stop = hasAllPureComponents(a,n) &&  (dmax <= d);

%[size(M,2) counter]
end
