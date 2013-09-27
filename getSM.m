function [M d0] = getSM(polysys,d,varargin)
% [SM d0] = getSM(polysys,d)
% --------------------------
% 
% This function makes the symbolic M matrix of given set of polynomial equations
% polysys and maximal total degree of d. The maximum total degree of the
% original set of equations d0 is also returned.
%
% M         =   Symbolic M matrix, cell of 
%
% d0        =   scalar, maximum total degree of the original set of polynomial
%               equations
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d         =   scalar, desired maximum total degree of matrix M
%
% CALLS
% -----
%
% getMon.m
%
% Kim Batselier, 2011-04-05

% number of variables
n = size(polysys{1,2},2);
% number of equations
n_eq = size(polysys,1);
% this vector will contain the degree of each equation
dorig = zeros(n_eq,1);

% get the full base of monomials for the original set of equations, need
% maximum degree present in the set of equations
d0 =0;

r= n_eq;    % initialize number of rows of M
for i = 1 : n_eq    
    dorig(i,1) = max(sum(polysys{i,2},2));
    if max(sum(polysys{i,2},2)) > d0
        d0 = max(sum(polysys{i,2},2));
    end
    r = r + nchoosek(d-dorig(i)+n,n)-1;     % update number of rows
end

% check given degree argument
if d < d0
    warning('You have provided a degree smaller than the maximum total degree of the given set of equations. Setting the degree to the maximum total degree in the set of equations...')
    d = d0;
end

% first allocate memory for the M matrix to speed up things
M = cell(r,1);
rowcounter = 1; 

%determine up to which degree we can do all equations together
dshared = min(d-dorig);

% do shared part

if dshared >= 0
    
    addBase = getMon(dshared,n);
    
    for j = 1 : size(addBase,1)     % for each shift
        
        for i =1: n_eq    % for each equation
                      
            M{rowcounter} = [exp2str(addBase(j,:)) ' f' num2str(i)];
            rowcounter = rowcounter + 1;
        end
    end    
end

% determine which degrees left
dleft = (d-dorig)-dshared;

% now do remaining parts

for i =1: n_eq    % for each equation
    
    % determine the monomials we additionally need to multiply with
    addBase = getMon(dshared+dleft(i),n,dshared+1);
    
    for j = 1 : size(addBase,1)     % for each shift
         
        M{rowcounter} = [exp2str(addBase(j,:)) ' f' num2str(i)];
        rowcounter = rowcounter + 1;
    end
end

    function str = exp2str(exp)
        
        str = [];
        for k = 1 : length(exp)
            if exp(k) == 0
            elseif exp(k) == 1
                str = [str ' x' num2str(k)];
            else
                str = [str ' x' num2str(k) '^' num2str(exp(k)) ];
            end
        end
    end
end
