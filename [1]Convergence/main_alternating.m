clc;clear all;
%cvx_solver sedumi
cvx_solver sdpt3
%cvx_solver mosek
%cvx_solver scs
%cvx_solver_settings('MAX_ITERS', 10^4, 'SCALE', 100);
%cvx_quiet(true)

%% Problem Data (I)
channel_generating=0;
save_data=1;
LC=1; %# loops for channel realizarions

CBF=true;
GSBF_Heuristic=true;
GSBF_Alternating=false;
Exhaustive=true;

%%
Area=10^4;
L=10; %Number of RAUs
N1=2;  % Antennas of each RAU

M=3; %Number of Multicast Groups
K=2; %Number of Mobile Users in Each Multicast Group

%Q=[6:-2:0]';  %QoS in dB
Q=4;
epsilon=0.05;  %Shape of the errors

%%
N_set=N1*ones(L,1); % Set of Antennas for all the RAU
K_set=K*ones(M,1); %Set of Numbers of Mobile Users in the Multicast Groups

Pc=5.6*ones(L,1);  %power consumption of the fronthaul link
amcoeff=1/4; %amplifier efficiency coefficient


%% Problem Data for Params (II)
params.K_set=K_set;
params.delta_set=ones(M,1);

params.N_set=N_set;
params.P_set=10^(0)*ones(L,1);   %set of transmit power constraints for all the RAUs

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
TotalPower_CBF_temp=0; TotalPower_GSBFHeuristic_temp=0;...
TotalPower_GSBFAlternating_temp=0; TotalPower_Exhaustive_temp=0; %network power consumption

TransmitPower_CBF_temp=0;  TransmitPower_GSBFHeuristic_temp=0;... 
TransmitPower_GSBFAlternating_temp=0;  TransmitPower_Exhaustive_temp=0; %Total Transmit Power consumption

A_number_Heuristic_temp=0;  A_number_Alternating_temp=0;...
A_number_Exhaustive_temp=0;  %Optimal number of active RAUs

SU_counter=0;

if channel_generating==true
for ss=1:LC   %generate channel
    for m=1:M
        [H(:,m,:,ss), Theta(:,:,m,:,ss)]=channel_realization(L, K_set(m), N_set, Area, epsilon);  %NxMxK: Estimated Channel
    end
end
save('H.mat','H');
save('Theta.mat','Theta');
end

 load('H.mat');
 load('Theta.mat');

lp=1;
params.H=H(:,:,:,lp);  %NxMxK channel matrix
params.Theta=Theta(:,:,:,:,lp);  %NxNxMxK channel covariance matrix
params.r_set=10^(Q/10)*ones(M,1);  %set of SINR thresholds  for each group 
    

mu=1/L.*ones(L,1);
TT=21;
value=zeros(TT,1);

for tt=1:TT

%% Compute Q: n+1 iteration
params.Inactive_index=[]; 
Active_index=[1:L];
params.Active_index=[1:L];
params.weight=1./mu;
params.rankone=false; %return general rank solution
[Q_GSBF, feasible_GSBF] = powermin_cvx(params);

value_temp=zeros(L,1);

for l=1:L 
    temp1=0;
    temp1_fro=0;
    for m=1:M
        temp1=temp1+trace(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
        if tt==1
        temp1_fro=temp1_fro+norm(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m),'fro');
        else
        temp1_fro=temp1_fro+norm(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m)-Q_old(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m),'fro');
        end
    end
    value_temp(l)=temp1+10^(-3);
    value_temp_fro(l)=temp1_fro;
end

Q_old=Q_GSBF;

%% Set mu: n+1 iteration
mu=sqrt(value_temp)./sum(sqrt(value_temp));
%% Objective value
value(tt)=(1./mu)'*value_temp;
%value_power(:,tt)=value_temp;
value_power(tt)=sum(value_temp_fro);
end

%plot([0:TT-1],value_power,'b-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
plot([0:TT-1],value*(5.6*4*4),'b-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
hold on;
xlabel('Iteration','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Objective Value','fontsize',14,'fontweight','b','fontname','helvetica');

if save_data==true
    Convergence=zeros(TT,2);
    Convergence(:,1)=[0:TT-1];
    Convergence(:,2)=value*(5.6*4*4);
    save('Convergence.dat','Convergence','-ascii');
end
