function [Vsolution,feasible] = powermin_cvx(params)

%prob_to_socp: maps PARAMS into a struct of SOCP matrices
%input struct 'parms' has the following fields:
%params.L;   %'L': # RRHs
%params.K;    %'K': # MUs
%params.N_set;  %set of antennas at all the RRHs
%params.Inactive_index;  %Index of the inactive RRHs
%params.Active_index;   %index of the active RRHs
%params.delta_set; %set of noise covariance

%%%%%%%%%%%%%%Problem Instances%%%%%%%%%%%%%
%params.r_set;  %set of SINR thresholds
%params.H;  %Channel Realization
%params.P_set;   %set of transmit power constraints at all the RRHs

%%%%%%%%Problem Data%%%%%%%
K_set=params.K_set;   %Mx1 vector: Set of Numbers of Mobile Users in M Multicast Group
r_set=params.r_set;     %Mx1 vector: QoS Requirements of Mobile Users for each Muliticast Group: assuming all the users have the same QoS requirments in the same group
delta_set=params.delta_set;  %Mx1 vector: noise covariance: assuming same value in the same group

N_set=params.N_set;   %Lx1 vector: RAU antennas set
P_set=params.P_set;  %Lx1 vector: RAU transmit power set

Inactive_index=params.Inactive_index;  %Index of the inactive RRHs
Active_index=params.Active_index;       %Index of the inactive RRHs

H=params.H;  %NxMxK channel matrix: N=sum(N_set), M=length(K_set), K equals the number of mobile users in each group (assuming each group has the same number of mobile users)
Theta=params.Theta; %NxNxMxK channel covariance matrix

weight=params.weight;  %length(Active_index)x1: weigth for each group beamformer
rankone=params.rankone;  %reture rankone solution or not

%%%%%%%%CVX+SCS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin

variable Q(sum(params.N_set), sum(params.N_set), length(K_set)) hermitian semidefinite;   %Variable for N x N x M beamforming matrix for M groups
variable t(length(Active_index),1);  %slack variables for each beamforming matrix for the RAUs
variable lambda(length(K_set), K_set(1));  % MxK variables for the S-Lemma in the QoS constraints
minimize (weight'*t)  % group sparsity inducing minimization
subject to

%% RAUs Transmit Power Constraints
%% Active RAUs
for l=1:length(Active_index)  
    temp1=0;
    for m=1:length(K_set)
        temp1=temp1+trace(Q(sum(N_set(1:Active_index(l)-1))+1:sum(N_set(1:Active_index(l))),sum(N_set(1:Active_index(l)-1))+1:sum(N_set(1:Active_index(l))),m));
    end
    temp1<=P_set(Active_index(l));
    temp1<=t(l);
    t(l)>=0;
end

%% Inactive RAUs
for l=1:length(Inactive_index)  
    temp1=0;
    for m=1:length(K_set)
        temp1=temp1+trace(Q(sum(N_set(1:Inactive_index(l)-1))+1:sum(N_set(1:Inactive_index(l))),sum(N_set(1:Inactive_index(l)-1))+1:sum(N_set(1:Inactive_index(l))),m));
    end
    temp1==0;
end


%% Semidefinite Constraints
% for m=1:length(K_set)
%     Q(:,:,m)==semidefinite(sum(N_set));
% end

%% QoS Constraints
Q_sum=0;
for m=1:length(K_set)
    Q_sum=Q_sum+Q(:,:,m);
end

for m=1:length(K_set)
    for k=1:K_set(m)
        G=(1+r_set(m))*Q(:,:,m)-r_set(m)*Q_sum;
        
        [G, G*H(:,m,k); H(:,m,k)'*G, H(:,m,k)'*G*H(:,m,k)-r_set(m)*delta_set(m)^2]+...
        lambda(m,k).*[Theta(:,:,m,k), zeros(sum(N_set),1); zeros(1,sum(N_set)), -1]==hermitian_semidefinite(sum(N_set)+1);
    
        lambda(m,k)>=0;
    end
end
cvx_end

 
%% build output
%%
     if  strfind(cvx_status,'Solved') 
         feasible=true;
         
         if rankone==true
         Vsolution=Rankone(Q);
         else
         Vsolution=Q;
         end
         
     else
         feasible=false;
         Vsolution=[];
     end