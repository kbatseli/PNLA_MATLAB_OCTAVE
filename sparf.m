function [root varargout] = sparf(polysys)
% [root d c ns check cr digits] = sparf(polysys)
% ----------------------------------------------
% Sparse affine root-finding algorithm. Try the quick & dirty algorithm first.
%
% root      =   matrix, each row corresponds with a component, each column
%               with an affine solution
%
% d         =   scalar, optional, degree for which all solutions were found
%
% c         =   vector, optional, degrees of each pure component in reduced
%               monomial system
%
% ns        =   vector, contains indices of rows of kernel of the normal
%               set
%
% check     =   scalar, norm(M*K,1), to check whether found solutions
%               evaluate to zero
%
% cr        =   vector, contains estimates of right corank up to d
%
% digits    =   vector, provides an upper bound on the maximum number of
%               decimal digits which are correct. If negative then no
%               decimal digits are correct.
%
% CALLS
% -----
% 
% getM.m, getD0.m, reduceAsys.m, hasAllPureComponents.m, getCra.m,
% getMex.m, updateN.m
%
% REFERENCES
% ----------
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

n=size(polysys{1,2},2);
d=getD0(polysys);
M=getM(polysys,d,1);
[Q, R, P]=qr(M');
r=nnz(diag(R));
c(d)=size(M,2)-r;
N=Q(:,r+1:end);

clear Q R P
h=0;

while ~h
    
    % reduced canonical decomposition on basis for null space N
    checki = [1:size(N,1)]; %contains indices of monomials that need to be checked for linear independence
    a=[];
    br=[];
    tol=sum(size(N))*eps;
    counter=1;
    
    while counter <= length(checki)
        [Y, Sin Z]=svd(full(N([br checki(counter)],:)'));
        if size(Sin,2)==1
            sin=Sin(1,1);
        else
            sin=diag(Sin);
        end
        rs=sum(sin > tol);
        
        if (asin(Sin(min(size(Sin)),min(size(Sin)))) < tol) || (rs < size(Sin,2))
            a = [a checki(counter)];
		% remove all monomial multiples from checki(counter)
		di=sum(fite(checki(counter),n));
		if di<d
            multiplei=2:nchoosek(d-di+n,n); % indices of monomial multiples
            for j=1:length(multiplei)
                checki(checki==feti(fite(checki(counter),n)+fite(multiplei(j),n)))=[];
            end
        end
        else
           br=[br checki(counter)];
        end        
        counter = counter + 1;
    end           
    
    % reduce the canonical decomposition and check pure powers
    Asys=cell(length(a),2);
    for i=1:length(a)
        Asys{i,2}=fite(a(i),n);
    end
    [h components] = hasAllPureComponents(Asys);
        
    
    if ~h ||  sum(components)-n <=0
        d=d+1;
        N=updateN(N,getMex(polysys,d,d-1,1),1);
        c(d)=size(N,2);
    else
        b=br;
        if c(d)>length(b) % check whether there are roots at infinity
            
            % determine index in N of 1st root at infinity            
            for i=a(end)+1:size(N,1)
                [Y, Sin Z]=svd(full(N([b i],:)'));
                if size(Sin,2)==1
                    sin=Sin(1,1);
                else
                    sin=diag(Sin);
                end
                rs=sum(sin > tol);
                
                if (asin(Sin(min(size(Sin)),min(size(Sin)))) < tol) || (rs < size(Sin,2))
                else
                    b=[b i];
                    break
                end
            end
        end
        
        % check whether shifting br takes us onto rows corresponding with
        % roots at infinity
        if length(br) > length(b) || (length(br)==length(b) && size(N,1)<=nchoosek(sum(fite(br(end),n))+n,n) ) || (length(b)>length(br) && sum(fite(b(length(br)+1),n)) == sum(fite(b(length(br)),n))+1)
            % shift took us to a root at infinity, need 1 more iteration
            h=0;
            d=d+1;
            N=updateN(N,getMex(polysys,d,d-1,1),1);
            c(d)=size(N,2);
        else
            % everything OK, proceed to construct the eigenvalue problem
            db   = nchoosek(sum(fite(br(end),n))+n,n);
            dbp1 = nchoosek(sum(fite(br(end),n))+n+1,n);
            [y,s, W]=svd(full(N(1:db,:)));
            Z=N(1:dbp1,:)*W;
            B=Z(1:db,1:length(br));
            
            A=zeros(size(B,1),size(B,2));
            % better to shift with random linear combination of components in
            % order to tackle multiplicities of a component
            coef = randn(1,n);  % take random components
            for j = 1 : n
                for i = 1:db            
                    A(i,:) = A(i,:) + coef(j)*Z(feti([zeros(1,j-1) 1 zeros(1,n-j)]+fite(i,n)),1:length(br));
                end
            end
            
            [V D]=eig(pinv(B)*A);                        
            
            K=Z(1:n+1,1:length(br))*V;
            K=K*diag(1./K(1,:));
            root=K(2:end,:);
            
            % check whether we found the roots
            K = makeRoot(getD0(polysys),root.');        
            for i = 1 : size(K,2)
                check(i) = norm(polysys2vec(polysys,getD0(polysys))*K(:,i)); % 
            end
            getal=0;for i=1:getD0(polysys),getal=getal+i*nchoosek(i+n-1,n-1);end        
            [check I] = sort(check);
            digits = floor(-log10(check/getal));
            root = root(:,I).';
            
            varargout{1}=d;
            varargout{2}=components;
            varargout{3}=br;
            varargout{4}=check;
            varargout{5}=c;
            varargout{6}=digits;
            
        end
        
        
    end
end

