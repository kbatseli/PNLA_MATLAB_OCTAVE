function [N nM tol]=sparseN(polysys,d)

s=size(polysys,1);
n=size(polysys{1,2},2);

for i=1:s
    polysys{i,1}=polysys{i,1}/norm(polysys{i,1},1);
end

% if M is sparse column compressed than accessing columns
% is faster than rows
M=getM(polysys,d,1); 
[p q]=size(M);

tol=max(p,q)*eps(sqrt(s));
nM=zeros(1,p);
nM(1)=norm(M(1,:));

nindex=find(M(1,:));
zindex=setdiff(1:q,nindex);

% initial way of constructing dense N
% N=zeros(q,q-1);
% for i=1:length(zindex)
%     N(zindex(i),i)=1;
% end

% initalize sparse N, first time we need to make it explicit
N=sparse(zindex,1:length(zindex),ones(1,length(zindex)),q,q-1);

pivotI=nindex(find(abs(M(1,nindex))==max(abs(M(1,nindex))),1));
M(1,:)=M(1,:)/M(1,pivotI);
tempI=setdiff(nindex,pivotI);

for j=1:length(tempI)
    N(pivotI,length(zindex)+j)=M(1,tempI(j));
    N(tempI(j),length(zindex)+j)=-M(1,pivotI);       
end
M=M*N;
[p q]=size(M);

for i=2:p
    nM(i)=norm(M(i,:));
    if nM(i) > tol
        nindex=find(M(i,:));
        zindex=setdiff(1:q,nindex);
        tempN=zeros(q,q-1);
        for j=1:length(zindex)
            tempN(zindex(j),j)=1;
        end
        %     temp=full(sparse(zindex,1:length(zindex),ones(1,length(zindex)),q,length(zindex)));
        pivotI=nindex(find(abs(M(i,nindex))==max(abs(M(i,nindex))),1));
        tempI=setdiff(nindex,pivotI);
        M(i,:)=M(i,:)/M(i,pivotI);
        for j=1:length(tempI)
            tempN(pivotI,length(zindex)+j)=M(i,tempI(j));
            tempN(tempI(j),length(zindex)+j)=-M(i,pivotI);       
        end
        M=M*tempN;
        [p q]=size(M);
        N=N*tempN;    
    end
end