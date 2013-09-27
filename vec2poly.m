function poly = vec2poly(vec,n,varargin)
% poly = vec2poly(vec,n,maple)
% ----------------------------
% Converts a row vector into a string representing the polynomial.
%
% poly      =   cell, each row contains the polynomial corresponding with
%               the rows of vec
%
% vec       =   matrix, each row is interpreted as a coefficient vector of
%               a polynomial
%
% n         =   scalar, number of indeterminates
%
% maplestr  =   boolean, if set to 1 will output a string which can be used
%               in maple, default=0.
%
% CALLS
% -----
% 
%
% Kim Batselier, 2011-08-01

if nargin == 2
    maplestr = 0;
else
    maplestr = varargin{1};
end
[p q] = size(vec);

poly = cell(p,1);

for i = 1 : p
    coefI = find(vec(i,:));
    ncoef = length(coefI);
    for j = 1 : ncoef
        if vec(i,coefI(j)) < 0
            signstr = ' ';
        else
            signstr = ' + ';
        end
        if ~maplestr
            poly{i,1} = [poly{i,1} signstr num2str(vec(i,coefI(j))) ' ' exp2str(frte(n,coefI(j))) ];
        else
            poly{i,1} = [poly{i,1} signstr num2str(vec(i,coefI(j))) '*' exp2str(frte(n,coefI(j))) ];
        end
    end
end


   function str = exp2str(exp)
        
        str = [];
        for z = 1 : length(exp)
            if exp(z) == 0
            elseif exp(z) == 1
                str = [str ' x_' num2str(z)];
            else
                str = [str ' x_' num2str(z) '^' num2str(exp(z)) ];
            end
        end
    end

end
        



