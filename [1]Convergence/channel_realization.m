function [H, Theta] =channel_realization(L, K, N_set, Area, epsilon)
% Generate channel matrix

%INPUT:
%L      =# RRHs
%K      =# MUs
%N_set    =# antennas at each RRH
%Area      =length of the area
%epsilon  = spherical error model

%OUTPUT:
% H     = (sum(N_set) x K) channel matrix
% Theta   =  (sum(N_set) x sum(N_set) x K) channel matrix

%%%%%%%%%%%%%%%%Network Realization%%%%%%%%%%%%%%%%%%%%%%%%
B_position=Area.*(rand(2,L)-0.5);  %%RRH positions
U_position=Area.*(rand(2,K)-0.5);  %% user positions

% B_position(1,1)=0;  B_position(2,1)=0;
% B_position(1,2)=Area/4;  B_position(2,2)=Area/4;
% B_position(1,3)=Area/4;  B_position(2,3)=-Area/4;
% B_position(1,4)=-Area/4;  B_position(2,4)=Area/4;
% B_position(1,5)=-Area/4;  B_position(2,5)=-Area/4;

 
%  U_position=zeros(2,K);
%  for k=1:K
%      if k==1
%          U_position(:,1)=Area.*(rand(2,1)-0.5); 
%      else
%          U_position(:,k)=U_position(2,1)+10*k;
%      end
%  end
 
 %%%%%Generate Large-Scale Fading%%%%%%%%%%%%%
for k=1:K
    for l=1:L
                   d=(norm(B_position(:,l)-U_position(:,k))+10);
                   %D(l,k)=5.6234*10^(5.355)*d^(-1.88)*sqrt(exp(normrnd(0,6.3)));
                   D(l,k)=1;
    end
end

%%%%%%Generate Small-Scale Fading%%%%%%%%%%%%%
for k=1:K
    for l=1:L
        %% Geneate estimated channel
        H(sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),k)=D(l,k)*(normrnd(0,1/sqrt(2),N_set(l),1)+i*normrnd(0,1/sqrt(2),N_set(l),1));   %noise normalized to 1
        %% Generate estimation errors
        %Theta(sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),k)=(epsilon^(-2)*D(l,k)^(-2))*(eye(N_set(l))); 
        Theta(sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),k)=(epsilon^(-2))*(eye(N_set(l))); 
    end
end