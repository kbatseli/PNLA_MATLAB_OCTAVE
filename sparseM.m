function M=sparseM(polysys,d,varargin)

dorig=getDorig(polysys);
neq= size(polysys,1);
n=size(polysys{1,2},2);
ncoef=zeros(neq,1);     % number of nonzero entries of each polynomial
nshifts=zeros(1,neq);   % number of rows in M(d) corresponding with each polynomial

if isempty(varargin)
 % no optional argument, default Macaulay matrix
 shiftindices = cell(neq,1);
 for i=1:neq
  ncoef(i) = length(polysys{i,1});
  shiftindices{i}=[1:nchoosek(d-dorig(i)+n,n)];
  nshifts(i)=length(shiftindices{i});
 end
else
 % optional argument specifies up to which degree we already have the Macaulay matrix
 dmin=varargin{1};
 % check optional input
 if dmin > d
     error('The extra argument cannot be smaller than d')
 end
 
 shiftindices = cell(neq,1);
 for i=1:neq
  ncoef(i) = length(polysys{i,1});
  shiftindices{i} = [nchoosek(dmin-dorig(i)+n,n)+1 : nchoosek(d-dorig(i) +n,n)];
  nshifts(i)=length(shiftindices{i});
 end
end

% specify the Macaulay matrix in triplet form
% total number of nonzero elements
nnz=nshifts*ncoef;
I=zeros(1,nnz);
J=zeros(1,nnz);
V=zeros(1,nnz);

% determine dimensions of Macaulay matrix
p=sum(nshifts);
q=nchoosek(d+n,n);

IVcounter=1;
Jcounter=1;
for i=1:neq
    for j=1:length(shiftindices{i})
        % coefficients in V
        V(IVcounter:IVcounter+ncoef(i)-1)=polysys{i,1};   
        % row indices in I
        I(IVcounter:IVcounter+ncoef(i)-1)=feti(polysys{i,2}+ones(ncoef(i),1)*fite(shiftindices{i}(j),n));
        % column indices in J
        J(IVcounter:IVcounter+ncoef(i)-1)=Jcounter*ones(1,ncoef(i));
        Jcounter=Jcounter+1;
        IVcounter=IVcounter+ncoef(i);
 end
end 

M=sparse(J,I,V,p,q);
end
