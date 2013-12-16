function [M nzmax] = getM(polysys,d,varargin)
% [M nzmax]  = getM(polysys,d,sparseM|factorvec)
% ----------------------------------------------
% 
% This function makes the M matrix of given set of polynomial equations
% polysys and maximal total degree of d. The maximum total degree of the
% original set of equations d0 is also returned. When the given degree d is
% smaller than d0 then the matrix of polynomials with degree <= d will be
% made. If there are no such polynonials then an empty matrix is returned.
%
% M         =   Matrix M of degree d corresponding with the set of equations
%               specified in polysys
%
% nzmax     =   scalar, total number of nonzero elements in M(d)
%
% polysys   =   cell containing coefficients and monomials exponents of the
%               set of polynomial equations.
%
% d         =   scalar, desired maximum total degree of matrix M
%
% sparseM   =   boolean, optional: set to 1 if M needs to be sparse,
%               default = 0
%
% factorvec =   cell, optional: polysys object of a factor. Only to be used
%               for the case M(d) is used for multiplication
%
% CALLS
% -----
%
% getMon.m, feti.m, getD0.m
%
% Kim Batselier, 2009-11, updated: removed one for loop to determine number
% of rows of M

% number of variables
n = size(polysys{1,2},2);
% number of equations
n_eq = size(polysys,1);
% this vector will contain the degree of each equation
dorig = zeros(n_eq,1);

% get the full base of monomials for the original set of equations, need
% maximum degree present in the set of equations
d0 =0;

if nargin < 3
    sparseM = 0;
    factorvec=[];
elseif nargin < 4 
    sparseM = varargin{1};
    factorvec=[];
else
    sparseM = varargin{1};
    if ~isempty(varargin{2})
        % the factor is ONLY for when M is used to do multiplications,
        % then we do not need to consider zero entries of the factor
        if iscell(varargin{2})
            % convert to vector
            factorvec = getM(varargin{2},getD0(varargin{2}),1);
            dfactor = getD0(varargin{2});
        else
            factorvec = varargin{2};
            dfactor = deg(factorvec,n);
        end
        % store nonzero entries in factor
        rowsToDo = find(factorvec);
    end
end

ncoef = zeros(1,n_eq);
for i = 1 : n_eq    
    ncoef(i) = length(polysys{i,1});
    dorig(i,1) = max(sum(polysys{i,2},2));
    if max(sum(polysys{i,2},2)) > d0
        d0 = max(sum(polysys{i,2},2));
    end    
end

% check given degree argument
if d < d0
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
        M = [];
        nzmax = 0;
        return
    end
end

%determine up to which degree we can do all equations together
dshared = min(d-dorig);

% determine which degrees cannot be done together
dleft = (d-dorig)-dshared;

r= n_eq;    % initialize number of rows of M
nzmax = 0;
if isempty(factorvec)
    for i = 1 : n_eq
        r = r + nchoosek(d-dorig(i)+n,n)-1;     % update number of rows
        nzmax = nzmax + ncoef(i)*nchoosek(d-dorig(i)+n,n);
    end
else
    r = nchoosek(dfactor+n,n);
    nzmax = length(rowsToDo)*ncoef(1); % count only rows we fill in
end

% number of columns of M
c = nchoosek(d+n,n);

% base of monomials V corresponding with original set of equations
% maxBase = getMon(d,n,ordering);

% first allocate memory for the M matrix to speed things up
if sparseM
    M = sparse([],[],[],r,c,nzmax); 
else
    M = zeros(r,c);    
end

rowcounter = 1; 

if ~isempty(factorvec)
  % we know at this point that dshared > 0 and dleft=0 
  
  % only do the rows that need to be done
  for i = 1 : length(rowsToDo)
      
      col = zeros(1,size(polysys{1,2},1));
      
      for j = 1 : size(polysys{1,2},1) % for each monomial in the equation
         col(j) = feti(fite(rowsToDo(i))+polysys{1,2}(j,:));
      end
      
      M(rowsToDo(i),col) = polysys{1,1};
  end
else
    % do shared part
    if dshared >= 0
    
    addBase = getMon(dshared,n);
    
    for j = 1 : size(addBase,1)     % for each shift
        
        for i =1: n_eq    % for each equation
            col = zeros(1,size(polysys{i,2},1));
            for k = 1 : size(polysys{i,2},1) % for each monomial in the equation
                
                col(k) = feti(addBase(j,:)+polysys{i,2}(k,:));            
            end
            
            M(rowcounter,col) = polysys{i,1};
            rowcounter = rowcounter + 1;
        end
    end    
    end
end



% now do remaining parts

for i =1: n_eq    % for each equation
    
    % determine the monomials we additionally need to multiply with
    addBase = getMon(dshared+dleft(i),n,dshared+1);
    
    for j = 1 : size(addBase,1)     % for each shift
        col = zeros(1,size(polysys{i,2},1));
        for k = 1 : size(polysys{i,2},1) % for each monomial in the equation
            
            col(k) = feti(addBase(j,:)+polysys{i,2}(k,:));            
        end        
        
        M(rowcounter,col) = polysys{i,1};
        rowcounter = rowcounter + 1;
%         col = [];
    end
end



end