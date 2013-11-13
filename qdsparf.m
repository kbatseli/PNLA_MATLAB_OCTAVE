function [root varargout] = qdsparf(polysys)
% [root d c ns check cr digits] = qdsparf(polysys)
% ------------------------------------------------
% Quick & Dirty Sparse Affine Root-Finding algorithm. Uses SuiteSparseQR.
%
% root      =   matrix, each row corresponds with a component, each column
%               with an affine solution
%
% d         =   scalar, optional, degree for which all solutions were found
%
% c         =   vector, optional, degrees of each pure component in reduced
%               monomial system
%
% ns        =   vector, contains indices of rows of kernel of the normal
%               set
%
% check     =   scalar, norm(M*K,1), to check whether found solutions
%               evaluate to zero
%
% cr        =   vector, contains estimates of right corank up to d
%
% digits    =   vector, provides an upper bound on the maximum number of
%               decimal digits which are correct. If negative then no
%               decimal digits are correct.
%
% CALLS
% -----
% 
% getM.m, spqr.m, getD0.m, getDreg.m, reduceAsys.m, hasAllPureComponents.m,
% getCra.m, getMex.m
%
% REFERENCES
% ----------
%
% T.A. Davis. Algorithm 9xx, SuiteSparseQR: multifrontal multithreaded
% rank-revealing sparse QR factorization. submitted to ACM Transactions on
% Mathematical Software, V:1–24, 20xx
%
% Kim Batselier, 2011-11-14, update 2011-11-15: builds up M iteratively
% with getMex.m

root = [];
d = getD0(polysys);
n = size(polysys{1,2},2);
% dmax = getDreg(polysys)+10; % conform 1, what is a good upper bound? STLS Philippe, dmax+2!

for i=1:size(polysys,1)
    polysys{i,1}=polysys{i,1}/norm(polysys{i,1});
end

h = 0;

% construct initial M matrix
[Ms nztemp] = getM(polysys,d-1,1);

while ~h
    
    % extend M with extra rows and columns
    [p q] = getMDim(polysys,d);
    temp = Ms;
    [Mex nzmex] = getMex(polysys,d,d-1,1);
    Ms = spalloc( p,q ,nztemp+nzmex);
    Ms(1:size(temp,1),1:size(temp,2)) = temp;
    Ms(size(temp,1)+1:end,:) = Mex;
    clear Mex
       
    % first QR of M'
    [Qs, Rs,Ps] = spqr(Ms',struct('permutation','vector'));    
    % determine rank of M
    rankMs = length(find(diag(Rs)));    
    % get basis for kernel
    Vs = Qs(:,rankMs+1:end);    
    crs(d) = size(Vs,2);    
    clear Qs Rs
        
    % second QR of Vs'
    [Qv Rv Pv] = spqr(Vs',struct('Q','discard','permutation','vector'));

    % normal set for both affine and roots at infinity
    rows = sort(Pv(1:crs(d) ));
    
    % find corresponding linear independent columns of M
    acol = sort(Pv(crs(d)+1:end));    
    for i = 1 : length(acol)
        Asys{i,1} = 1;
        Asys{i,2} = fite(acol(i),size(polysys{1,2},2));
    end    
    % find reduced monomial system
    Asys = reduceAsys(Asys);
    
    [h c] = hasAllPureComponents(Asys);
    
    if sum(c)-n <0 % check whether we can compute an affine normal set
        ns =[];
    else
        ns = getCra(Asys,sum(c)-n); % normal set for only affine roots
    end
    
    if ~h
        d = d + 1;        
    else
        % normal set for affine roots
%         ns = getCra(Asys,sum(c)-n);
        
        if feti([1 zeros(1,n-1)]+fite(ns(end),n)) > size(Vs,1)
            h = 0;
            d = d+1;
        else
        
        % initialize root
        root = zeros(n,length(ns));
        check = zeros(1,length(ns));
        dns = sum(fite(ns,n));  % highest degree of element in affine normal set
        
        % third QR of Vs(ns,:) to determine which columns we'll need for
        % constructing the eigenvalue problem
        [Q R P] = spqr(Vs(ns,:),struct('Q','discard','permutation','vector'));
        
 
        % which columns of Vs do we need to make the eigenvalue problem
        columns = sort(P(1:length(ns)));
        
       B = full(Vs(ns,columns));
       % At this point B should be a square full rank matrix!
       if rank(B) ~= size(B,2)
            disp('warning: B matris is not of full rank')
            [size(B,2) rank(B)]
       end
       A = zeros(size(B,1),size(B,1));
        
        % better to shift with random linear combination of components in
        % order to tackle multiplicities of a component
       coef = randn(1,n);  % take random components
       for j = 1 : n
        for i = 1:size(B,1)            
%                 A(i,:) = A(i,:) + coef(j)*Vs(feti([zeros(1,j-1) 1 zeros(1,n-j)]+fite(i,n)),columns);
               A(i,:) = A(i,:) + coef(j)*Vs(feti([zeros(1,j-1) 1 zeros(1,n-j)]+fite(ns(i),n)),columns);
            end        
        end

       [Vb Db] = eig(A,B);

	if rank(Vb) ~= size(Vb,2)
		disp('warning: eigenvectors not of full rank')
	end
        
        % check whether we found the roots
        Ks = Vs(1:n+1,columns)*sparse(Vb);
        Ks = Ks * diag(1./Ks(1,:)); 
        root = Ks(2:2+n-1,:);
        K = makeRoot(getD0(polysys),root.');        
        for i = 1 : size(K,2)
            check(i) = norm(polysys2vec(polysys,getD0(polysys))*K(:,i)); % 
        end
        getal=0;for i=1:getD0(polysys),getal=getal+i*nchoosek(i+n-1,n-1);end        
        [check I] = sort(check);
        digits = floor(-log10(check/getal));
        root = root(:,I).';
        varargout{1}= d;
        varargout{2}= c;
        varargout{3} = ns;
        varargout{4} = check;
        varargout{5} = crs;
        varargout{6} = digits;
        
        end
    end    
end

if ~h
    error('Something went wrong. Please call the SMC hotline.')
end

end

