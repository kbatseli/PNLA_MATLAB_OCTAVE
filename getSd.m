function Sd = getSd(polysys,d)
% Sd = getSd(polysys,d)
% ---------------------
% 
% This function returns a matrix Sd such that the rows of Sd are only the
% shifts of the polynomial system polysys to degree d.
%
% Sd        =   matrix, rows correspond with only the shifts to degree d
%
% d         =   scalar, desired degree of matrix Sd
%
% CALLS
% -----
%
% getMon.m, feti.m

% number of variables
n = size(polysys{1,2},2);
% number of equations
n_eq = size(polysys,1);
% this vector will contain the degree of each equation
dorig = zeros(n_eq,1);

% determine degree of each equation
for i = 1 : n_eq
    dorig(i,1) = max(sum(polysys{i,2},2));
end

d0 = max(dorig);

% check given degree argument
if d < d0
    error(['You have provided a degree smaller than the minimal degree of M. Argument d should be at least ' num2str(d0)])
end

% initialize number of rows of M
r = 0;
for i = 1 : n_eq
    r = r + nchoosek(d-dorig(i)+n-1,n-1);
end

% number of columns of M
c = nchoosek(d+n,n);

% first allocate memory for the M matrix to speed up things
Sd = zeros(r,c);
rowcounter = 1; 

% determine degrees of shifts
dshift = d-dorig;

for i =1: n_eq    % for each equation
    
    % determine the monomials we additionally need to multiply with
    addBase = getMon(dshift(i),n,dshift(i));
    
    for j = 1 : size(addBase,1)     % for each shift
        
        for k = 1 : size(polysys{i,2},1) % for each monomial in the equation
            
            col(k) = feti(addBase(j,:)+polysys{i,2}(k,:));            
        end        
        
        Sd(rowcounter,col) = polysys{i,1};
        rowcounter = rowcounter + 1;
        col = [];
    end
end


end
