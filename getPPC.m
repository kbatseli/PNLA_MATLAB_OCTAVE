function ppc = getPPC(i,d,n)
% ppc = getPPC(i,d,n)
% -----------------
% Returns a coefficent vector corresponding with a polynomial in pure
% components of x_i of degree d in n unknowns.
%
% ppc       =   row vector, coefficent vector corresponding with a
%               polynomial in pure components of x_i of degree d in n
%               unknowns
%
% i         =   scalar, index of the unknown (eg x_1 or x_2)
%
% d         =   scalar, degree 
%
% n         =   scalar, total number of unknowns
%
% CALLS
% -----
% 
% getMon.m
%
% Kim Batselier, 2010-04-22

mon = getMon(d,n);

smon = sum(mon(:,[1:i-1 i+1:n]),2);

% index of pure components locations in monomial basis
ipc = find(smon==0);

% how many pure components are there
npc = length(ipc);

ppc = zeros(npc+1,size(mon,1));

for j = 1 : npc
   ppc(j,ipc(j)) = 1;
end

ppc(end,:) = (smon==0)';