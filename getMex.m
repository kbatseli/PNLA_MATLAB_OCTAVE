function [Mex nzmax]= getMex(polysys,d,dmin,varargin)
% Mex = getMex(polysys,d,dmin,sparseM)
% ------------------------------------
% 
% Returns the extra rows that need to be added to M(dmin) in order to
% construct M(d). Suitable for recursive updating of M.
%
% Mex       =   matrix, contains the extra rows to go from M(dmin) to 
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d         =   scalar, desired maximum total degree of matrix M
%
% dmin      =   scalar, degree of already built-up M matrix
%
% sparseM   =   boolean, optional: set to 1 if M needs to be sparse,
%               default = 0
%
% CALLS
% -----
%
%
% Kim Batselier, 2011-11-15

% check given degree argument
if d < dmin
    error('d should be larger than dmin')
end

if nargin < 4
    sparseM = 0;
else
    sparseM = varargin{1};
end

% number of variables
n = size(polysys{1,2},2);
% number of equations
n_eq = size(polysys,1);
% this vector will contain the degree of each equation
dorig = getDorig(polysys);

if d < max(dorig)
    % need to consider only the polynomials of degree d or smaller
    indices = find(dorig <= d);
    if ~isempty(indices)
        temp = cell(length(indices),2);
        for i = 1 : length(indices)
            temp{i,1} = polysys{indices(i),1};
            temp{i,2} = polysys{indices(i),2};
        end
        clear polysys
        polysys = temp;
        dorig = dorig(indices,1);
        n_eq = length(indices);
    else
        Mex = [];
        nzmax = 0;
        return
    end
end

r = 0; % initialize number of rows of Mex
nzmax = 0; % number of nonzero elements
ncoef = zeros(1,n_eq);
shifts = zeros(n_eq,d-dmin);
for i = 1 : n_eq    
    ncoef(i) = length(polysys{i,1});
%     dorig(i,1) = max(sum(polysys{i,2},2));
    shifts(i,:) = sort([d:-1:dmin+1]-dorig(i)); % get degree of shift monomials for each polynomial
    for j = 1 : length(shifts(i,:))
        r = r + nchoosek(shifts(i,j)+n-1,n-1); 
        nzmax = nzmax + ncoef(i)*nchoosek(shifts(i,j)+n-1,n-1); 
    end
end

% number of columns of M
c = nchoosek(d+n,n);

if n_eq > 1
    % find intersection of shifts
    shared = intersect(shifts(1,:),shifts(2,:));
    for i = 3:size(shifts,1)
        shared = intersect(shared,shifts(i,:));
    end
    % get remaining separate shifts
    for i = 1 : n_eq
        if ~isempty(setdiff(shifts(i,:),shared))
            separate(i,:) = setdiff(shifts(i,:),shared);
        else
            separate = [];
        end
            
    end
else
    % only 1 polynonial
    shared = shifts;
    separate = [];
end

% first allocate memory for the M matrix to speed up things
if sparseM
    Mex = sparse([],[],[],r,c,nzmax); 
else
    Mex = zeros(r,c);    
end

rowcounter = 1; 
% do shared part
if ~isempty(shared)
    
    addBase = getMon(max(shared),n,min(shared));
    
    for j = 1 : size(addBase,1)     % for each shift
        
        for i =1: n_eq    % for each equation
            col = zeros(1,size(polysys{i,2},1));
            for k = 1 : size(polysys{i,2},1) % for each monomial in the equation
                
                col(k) = feti(addBase(j,:)+polysys{i,2}(k,:));            
            end
            
            Mex(rowcounter,col) = polysys{i,1};
            rowcounter = rowcounter + 1;
%             col = [];
        end
    end    
end

if ~isempty(separate)
    for i = 1 : n_eq
        
        addBase = getMon(max(separate(i,:)),n,min(separate(i,:)));
        
        for j = 1 : size(addBase,1)
            col = zeros(1,size(polysys{i,2},1));
            for k = 1 : size(polysys{i,2},1) % for each monomial in the equation
                col(k) = feti(addBase(j,:)+polysys{i,2}(k,:));            
            end

            Mex(rowcounter,col) = polysys{i,1};
            rowcounter = rowcounter + 1;
%             col = [];
            
        end
    end
end
