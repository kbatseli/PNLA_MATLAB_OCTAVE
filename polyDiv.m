function [q r a e rM] = polyDiv(psys,polysys)
% [q r a e rM] = polyDiv(psys,polysys)
% ------------------------------------
%
% Divides a multivariate polynomial p (psys) by the set of multivariate 
% polynomials polysys. Graded xel monomial ordering is always assumed. 
% This decomposes p into
%
% p = q + r
%
% where q lies in the row space of the Divisor matrix. The remainder r is
% unique since it is expressed in terms of the standard monomials.
%
% q         =   vector, contains coefficients of quotient polynomial
%
% r         =   vector, contains coefficients of remainder polynomial
%
% a         =   vector, expresses q as linear combination of rows of
%               Divisor matrix M
%
% e         =   vector, contains indices of linear indepent rows of
%               Divisor matrix M
%
% rM        =   vector, indicates which rows of the Divisor matrix are
%               needed in the linear combination to construct q = a*M(rM,:)
%
% p         =   cell, polysys cell for the multivariate polynomial p
%
% polysys   =   cell, contains coefficients and monomials exponents of the
%               set of polynomial equations
%
%
% CALLS
% -----
%
% getD.m, feti.m
%
% Kim Batselier, 2010-10-04, update 2011-11: uses sparse matrices now,
% 2012: removed the use of intRound.m
% 2013: uses default QR of Matlab, hence this routine doesn't work in Octave
% anymore
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
% 
%     You should have received a copy of the GNU Lesser General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

n = size(polysys{1,2},2);

np = size(psys{1,2},2);

if n ~= np
    error('Polynomial and polynomial system do not have the same amount of variables.')
end

% determine degrees of polynomial p and polynomial system polysys
d=max(sum(psys{1,2},2));
d0=max(sum(polysys{1,2},2));
for i=2:size(polysys,1)
    if max(sum(polysys{i,2},2)) > d0
        d0 =max(sum(polysys{i,2},2));
    end
end
    
% vectorize polynomial p
pindices=zeros(1,length(psys{1,1}));
for i=1:length(psys{1,1})
    pindices(i)=feti(psys{1,2}(i,:));
end
p=sparse(1,pindices,psys{1,1},1,nchoosek(d+n,n),length(psys{1,1}));

if d < d0
    q =[];
    r = p;
    a = [];
    e=[];
    rM=[];
    return
end

D = getD(polysys,p,1);    % construct divisor matrix

% limit the p vector to the number of columns of M
p=p(1:size(D,2));

% %% determine rank of Divisor matrix
[Q R P] = qr(D','vector');
% [Q R P] = spqr(D',struct('Q','matrix','permutation','vector')); %
% SuiteSparseQR
rankD = length(find(diag(R)));


% store indices of independent rows of Divisor matrix
e = sort(P(1:rankD));
% Remove linear dependent rows of D
D = D(e,:); 
V = Q(:,rankD+1:end);

clear Q R P

% Determine normal set that will span the remainder
[~, ~, Pv] = qr(V','vector');
% [Qv Rv Pv] = spqr(V',struct('Q','discard','permutation','vector'));
rowI = sort(Pv(1:size(V,2)));
B = sparse(1:size(V,2),rowI,ones(1,size(V,2)),size(V,2),size(D,2),size(V,2));

R = qr([B;D;p]');

% coefficients of D according to Q basis
RD = R(:,size(B,1)+1:end-1)';
% coefficients of p according to Q basis
Rp = R(:,end)';
clear R

RDW = RD(:,end-rankD+1:end);
RpW = Rp(:,end-rankD+1:end);

% inverse of RDW with backslash
[C R E]=qr(RDW,speye(size(RDW,1),size(RDW,1)));
RDWinv = E*(R\C);
[msg] = lastwarn;
if strcmp(msg,'Matrix is singular to working precision.')
    disp('Using pseudo-inverse')
    RDWinv = pinv(full(RDW));
    lastwarn('')
end

a = RpW*RDWinv;
q = a*D;

rM = find(a);
r = p-q;

end
