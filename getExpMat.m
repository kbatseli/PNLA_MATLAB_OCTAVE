function expmatrix = getExpMat(varargin)
% expmatrix = getExpMat(monsys) or getExpMat(a,n)
% -----------------------------------------------
% Returns a matrix where each row is an exponent vector of a given monomial
% system.
%
% expmatrix =   matrix, each row contains an exponent 
%
% monsys    =   cell, polysys representation of a monomial system
%
% a         =   vector, contains indices of reduced linear independent leading
%               monomials
%
% n         =   scalar, number of variables
%
% CALLS
% -----
%
% frte.m
%
% Kim Batselier, 2012-04, update 2012-09-18: accepts a,n as input now.

if iscell(varargin{1})
    monsys = varargin{1};
    n = size(monsys{1,2},2);
    n_mon = size(monsys,1);
    expmatrix = zeros(n_mon,n);
    % convert the exponent cells of syzsys into a matrix
    for i = 1 : n_mon
        expmatrix(i,:) = monsys{i,2};
    end
else
    a = varargin{1};
    n=varargin{2};
    expmatrix = zeros(length(a),n);
    for i=1:size(expmatrix,1)
       expmatrix(i,:)=fite(a(i),n);
    end
end
