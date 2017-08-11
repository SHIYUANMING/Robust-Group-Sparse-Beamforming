function [W]=Rankone(Q)
%INPUT:
% Q  =(L*N1) x (L*N1) x (K) matrix

%OUTPUT
%W   =L*N1 x K matrix

%Using Gaussian Randomization Technology

K=size(Q, 3);
N=size(Q,1);
W=zeros(N, K);

%% Gaussian Randomization
% for kk=1:K
%     %%%%%%check if Q(:,:,k) is a rank one matrix%%%%%%
%     [S, D]=eig(Q(:,:,kk));  % S: eigvector;   D %eigvalue
%     Index=find(diag(D)>10^(-4));
%     if length(Index)==1  %rank one
%         W(:,kk)=sqrt(D(Index,Index))*S(:,Index);
%     else
%         W(:,kk)=S*sqrtm(D)*randn(N,1);  %Gaussian Randomization
%     end   
% end

%% Maximum eigvalue
for kk=1:K
    %%%%%%check if Q(:,:,k) is a rank one matrix%%%%%%
    [S, D]=eig(Q(:,:,kk));  % S: eigvector;   D %eigvalue
    [value,Index]=sort(diag(D));
    W(:,kk)=sqrt(D(Index(length(diag(D))),Index(length(diag(D)))))*S(:,Index(length(diag(D)))); 
end
