function [Nnew tol] = updateN(N,M,varargin)
% [Nnew tol] = updateN(N,M,sparse)
% --------------------------------
% Updates the numerical kernel N of M(d) using the recursive
% orthogonalization theorem.
%
% Nnew          =   matrix, updated numerical kernel of M(d+1), number of
%                   columns is the dimension of the null space
%
% tol           =   scalar, default numerical tolerance used in the rank-test
%
% N             =   matrix, numerical basis for the kernel N of M(d)
%
% M             =   matrix, extra rows that need to be added to M(d) to
%                   construct M(d+1)
%
% CALLS
% -----
%
% Kim Batselier, 2013

switch nargin
    case 2
        sparseM=0;
    case 3
        if ~issparse(M)
            warning('Storage class of provided matrix was not sparse. Converting to sparse matrix....')
            M=sparse(M);
        end
        sparseM=1;
        defaultTol = 1;
    case 4
        if ~issparse(M)
            warning('Storage class of provided matrix was not sparse. Converting to sparse matrix....')
            M=sparse(M);
        end
        sparseM=1;
        defaultTol = 0;
    otherwise
        error('Minimum 2 input arguments are required');
end


% first determine the sizes of the different block matrices
%[q r] = size(U);
[q m] = size(N);
if q==0 || m==0
    q=0;
    m=0;
    N=[];
end
[dp q2] = size(M);

% construct the part from which we can determine the increase of rank
temp = [N'*M(:,1:q)';M(:,q+1:end)'];

if ~sparseM
    [Ut St Vt] = svd(temp);
    st = diag(St);
    tol = eps(st(1))*max(size(temp));
    dr = sum(st > tol);
    if dr < length(st) && st(dr)/st(dr+1) < 100
        warning(['Warning: approxi-rank gap is ' num2str(st(dr)/st(dr+1)) ', numerical rank might be ill-defined.']);
    end
else
    [Ut R P] = qr(temp);
    tol = 20*sum(size(temp))*eps;        
%     dr = nnz(diag(R));
    
    %% not sure to keep this in, heuristical rank-estimation
    diagR = sort(abs(diag(R)),'descend'); % do some kind of `singular value analysis', sort them from large to small
    diagR=diagR/diagR(1); % normalize such that first `singular value' is 1
    dr1 = nnz(diagR); % number of nonzero elements in diagR
    diagR = diagR(1:dr1); % remove the zero coefficiencts from diagR so that we don't divide by zero
    [rankgap dr2] = max(diagR(1:end-1)./diagR(2:end));
    if defaultTol
        dr = dr1;
    elseif dr1==size(temp,2) % temp is of full column rank, hence do not choose `maximal rank-gap'
        dr = dr1;
    else
        dr=dr2;  
    end
    
end

%Unew = [U N*Ut(1:m,1:dr);zeros(q2-q,r) Ut(m+1:end,1:dr)];
Nnew = [N*Ut(1:m,dr+1:end);Ut(m+1:end,dr+1:end)];

end
