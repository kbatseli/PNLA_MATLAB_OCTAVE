function SSd = getSSd(polysys,d)
% SSd = getSSd(polysys,d)
% -----------------------
% 
% This function returns a symbolic matrix Sd such that the rows of Sd are only the
% shifts of the polynomial system polysys to degree d.
%
% SSd        =   matrix, rows correspond with only the shifts to degree d
%
% d         =   scalar, desired degree of matrix Sd
%
% CALLS
% -----
%
% getMon.m

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

% first allocate memory for the M matrix to speed up things
SSd = cell(r,1);
rowCounter = 1; 

% determine degrees of shifts
dshift = d-dorig;

for i =1: n_eq    % for each equation
    
    % determine the monomials we additionally need to multiply with
    addBase = getMon(dshift(i),n,dshift(i));
    
        % for each monomial in the tempBase
        for k = 1 : size(addBase,1)
            
            SSd{rowCounter} = [exp2str(addBase(k,:)) ' f' num2str(i)];
            
            rowCounter = rowCounter +1;
        end

end

   function str = exp2str(exp)
        
        str = [];
        for z = 1 : length(exp)
            if exp(z) == 0
            elseif exp(z) == 1
                str = [str ' x' num2str(z)];
            else
                str = [str ' x' num2str(z) '^' num2str(exp(z)) ];
            end
        end
    end


end
