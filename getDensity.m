function [density mem] = getDensity(A,varargin)
% density = getDensity(polysys,d) OR density = getDensity(M)
% ----------------------------------------------------------
% Returns the density of the Macaulay matrix M(d). 
%
% density   =   scalar, density of Macaulay matrix M(d) in percent %
%
% mem       =   scalar, memory in bytes needed to store the full M(d)
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d			=	scalar, degree for which M(d) has to be constructed from polysys
%
% M 		=	matrix, Macaulay matrix (instead of providing polysys)
%
% CALLS
% -----
%
% getM.m
%
% Kim Batselier, 2012-08-03

density = 0;
mem=0;

if iscell(A)
% input argument is a polysys cell, need extra degree argument
	if isempty(varargin)
		error('Error: when providing a polysys cell you also need to provide a degree');
	elseif ~isscalar(varargin{1})
		error('Error: extra input argument should be a scalar');
    else
        d0=0;
        d=varargin{1}; % get degree
        n = size(A{1,2},2);
        n_eq = size(A,1); % get number of polynomials
        ncoef = zeros(1,n_eq); % contains for each polynomials number of nonzero coefficients
        for i = 1 : n_eq    
            ncoef(i) = length(A{i,1});
            dorig(i,1) = max(sum(A{i,2},2));
            if max(sum(A{i,2},2)) > d0
                d0 = max(sum(A{i,2},2));
            end
        end
        if d < d0
            % need to consider only the polynomials of degree d or smaller
            indices = find(dorig <= d);
            if ~isempty(indices)
                temp = cell(length(indices),2);
                for i = 1 : length(indices)
                    temp{i,1} = A{indices(i),1};
                    temp{i,2} = A{indices(i),2};
                end
                clear A
                A = temp;
                dorig = dorig(indices,1);
                n_eq = length(indices);
            else
%                 density = [];
                error(['Given polynomial system does not contain any polynomial of total degree ' num2str(d)])
%                 return
            end
        end
        r= n_eq; % contains total amount of rows of M(d)
        nzmax = 0; % contains total amount of nonzero coefficients
        for i = 1 : n_eq
            r = r + nchoosek(d-dorig(i)+n,n)-1;
            nzmax = nzmax + ncoef(i)*nchoosek(d-dorig(i)+n,n);
        end
% 		M = getM(A,varargin{1},1,1);
		density = 100 * nzmax/ (r*nchoosek(d+n,n));
        mem = (r*nchoosek(d+n,n))*8; % 64 bit = 8 bytes per value
	end
else
	density = 100 * nnz(A)/ numel(A);
    mem = numel(A)*8;
end

end
