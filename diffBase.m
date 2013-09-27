function y = diffBase(d,root,x)
% y = diffBase(d,root,x)
% ----------------------
%
% Applies the (partial) differential operator on a polynomial base vector
% of degree 'd' and evaluates it with 'root'. Graded xelicographic ordering
% is implicitly assumed. Multiple differentiation is also supported, x
% should then be a vector which indicates in which order there needs to be
% differentiated.
%
% y     =   column vector, contains the evaluated polynomial base vector
%
%           d^n1x1 ... d^nnxn  |
%           ------------------ | 
%           dx1^n1 ... dxn^nn  |x = root
%
% d     =   scalar, degree of the multivariate polynomial base vector
%
% root  =   row vector, used to evaluate differentiated base with
%
% x     =   row vector, index of the variable to which the differentation needs
%           to take place,
%
%                d        d        d
%           1 = ---, 2 = ---, 3 = ---, etc...
%               dx1      dx2      dx3  
%               
%           when a higher order differentiation is required then x is a row
%           vector indicating the order in which the variables need to be
%           differentiated, eg. x = [1 1 2] means first a 2nd order
%           differentiation to x1, then 1st order differentiation to x2.
%
% EXAMPLE
% -------
%
% second order derivative to y (= Dyy) of degree 6
% evaluated in the point (2,3):
%
% diffBase(6,[2 3],[2 2])
%
% first 1st order derivative to y, then 1st order derivative to x (= Dyx),
% degree 4 and evaluated in (-5,9):
%
% diffBase(4,[-5 9],[2 1])
%
% CALLS
% -----
%
% getMon.m
%
% Kim Batselier, 2010-01-28

n = size(root,2);
monBase = getMon(d,n);
l = length(monBase);
% coef = zeros(l,1);
Dn = size(x,2);

if size(x,2) ~= Dn
    error('The ''x'' argument provided should be a vector indicating for each differentiation step to which variable needs to be differentiated')
end

% check function inputs
if x > n
        error(['You cannot differentiate with respect to x' num2str(x) ', there are only ' num2str(n) ' variables.'])
end

coef = ones(l,1);

% derivative of the monomial base
DmonBase  = monBase;
for i = 1 : Dn
    if exist('indices','var')
        clear indices
    end
    % first run of the coefficients
    indices = find(DmonBase(:,x(i)) ~= 0);
    coef(indices) = coef(indices).*DmonBase(indices,x(i));

    DmonBase = DmonBase-[zeros(l,x(i)-1) ones(l,1) zeros(l,n-x(i))];
    
    % negative exponents are manually put to zero
    DmonBase = DmonBase.*(~(DmonBase(:,x(i)) < 0)*ones(1,n));   
end

temp = zeros(l,n);
for i = 1:length(indices)
    temp(indices(i),:) = root.^DmonBase(indices(i),:);
end

y = zeros(l,1);
y = coef.*prod(temp,2)./getDenom(x);

    function denom = getDenom(x)
        xmax = max(x);
        exp = zeros(1,xmax);
        for i = 1 : xmax
            exp(i) = length(find(x== i));            
        end
        denom = prod(factorial(exp));
    end
end