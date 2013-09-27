function vec = zeng2vec(zeng)
% vec = zeng2vec(zeng)
% --------------------
%
% Convert a zeng matrix into a coefficient vector.
%
% vec   =   vector, row vector containing coefficients
%
% zeng  =   matrix, each column corresponds with a term, first n rows are
%           exponent, last row is coefficient
%
% CALLS
% -----
%
% feti.m
%
% Kim Batselier, 2012-02-16

n = size(zeng,1)-1;
d = max(sum(zeng(1:n,:)));

vec = spalloc(1,nchoosek(d+n,n),size(zeng,2));

for i = 1 : size(zeng,2)
   vec(feti(zeng(1:n,i))) = zeng(end,i);
end

end
