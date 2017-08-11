clc;clear all;
%cvx_solver sedumi
cvx_solver sdpt3
%cvx_solver mosek
%cvx_solver scs
%cvx_solver_settings('MAX_ITERS', 10^4, 'SCALE', 100);
cvx_quiet(true)

%% Problem Data (I)
channel_generating=1;
LC=1; %# loops for channel realizarions

CBF=1;
GSBF_Alternating=1;
Baseline=1;

Exhaustive=0;
GSBF_Heuristic=0;


%%
Area=4*10^3;
L=8; %Number of RAUs
N1=2;  % Antennas of each RAU

M=5; %Number of Multicast Groups
K=2; %Number of Mobile Users in Each Multicast Group

Q=[8:-2:0]';  %QoS in dB
%Q=6;
epsilon=0.05;  %Shape of the errors

%%
N_set=N1*ones(L,1); % Set of Antennas for all the RAU
K_set=K*ones(M,1); %Set of Numbers of Mobile Users in the Multicast Groups

%Pc=5.6*ones(L,1);  %power consumption of the fronthaul link
Pc=5.6*ones(L,1)+[0:L-1]';  %power consumption of the fronthaul link
amcoeff=1/4*ones(L,1); %amplifier efficiency coefficient


%% Problem Data for Params (II)
params.K_set=K_set;
params.delta_set=ones(M,1);

params.N_set=N_set;
params.P_set=10^(0)*ones(L,1);   %set of transmit power constraints for all the RAUs

%%%%%%initialize results%%%%%%%%%%%%%%%%%%%%%%
TotalPower_CBF_temp=0; TotalPower_GSBFHeuristic_temp=0;TotalPower_GSBFAlternating_temp=0;...
TotalPower_Baseline_temp=0; TotalPower_Exhaustive_temp=0; %network power consumption

TransmitPower_CBF_temp=0;  TransmitPower_GSBFHeuristic_temp=0; TransmitPower_GSBFAlternating_temp=0;... 
 TransmitPower_Baseline_temp=0; TransmitPower_Exhaustive_temp=0; %Total Transmit Power consumption

A_number_Heuristic_temp=0;  A_number_Alternating_temp=0;...
A_number_Baseline_temp=0; A_number_Exhaustive_temp=0;  %Optimal number of active RAUs

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

for lp=1:LC
    
    SU_channel_temp=1;  %recode if the channel is feasible
    
    params.H=H(:,:,:,lp);  %NxMxK channel matrix
    params.Theta=Theta(:,:,:,:,lp);  %NxNxMxK channel covariance matrix
    
    for lq=1:length(Q)
        params.r_set=10^(Q(lq)/10)*ones(M,1);  %set of SINR thresholds  for each group 
        
 %% Coordinated Beamforming%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CBF==true

%% Problem Solving
params.Inactive_index=[]; 
params.Active_index=[1:L];
%params.weight=(1/amcoeff)*ones(length(params.Active_index),1);
params.weight=(1./amcoeff(params.Active_index));
params.rankone=true; %return rankone solution
[V_CBF, feasible_CBF] = powermin_cvx(params);

%% Output for the Particular QoS
if feasible_CBF==1
TotalPower_CBF(lq)=(1/amcoeff(1))*norm(V_CBF,'fro')^2+sum(Pc);
TransmitPower_CBF(lq)=(1/amcoeff(1))*norm(V_CBF,'fro')^2;
else
    SU_channel_temp=0;  % this channel if infeasible
    break;
end

end       
        

%% GSBF with Heuristic Group Sparsity Design%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if GSBF_Heuristic==true
    
    %% Compute Approximated Group Sparse Beamformer
    params.Inactive_index=[]; 
    params.Active_index=[1:L];
    %params.weight=(Pc(1)/amcoeff*length(params.Active_index))*ones(length(params.Active_index),1);
    params.weight=Pc(params.Active_index)./amcoeff(params.Active_index); 
    params.rankone=false; %return general rank solution
    [Q_GSBF, feasible_GSBF] = powermin_cvx(params);
    
    if feasible_GSBF==false
        SU_channel_temp=0;  %This channel is infeasible
        break;
    else  
        %% RRH Ordering Based on the Approximated Group Sparse Beamformer
        Value_GSBF=zeros(L,1); 
        for l=1:L 
            temp1=0;
            for m=1:M
                temp1=temp1+trace(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
            end
            Value_GSBF(l)=temp1;
        end
        
        %% Sparsity Parttern Based RRH Ordering 
        %[BS,BS_index]=sort(Value_GSBF);
        
        %% Proposed RRH Ordering
        ChannelGain=zeros(L,1);
        for l=1:L
            ChannelGain(l)=norm(vec(params.H(sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),:,:)),'fro')^2;
        end
        
        %[BS,BS_index]=sort(ChannelGain.*Value_GSBF);
        [BS,BS_index]=sort(((ChannelGain./Pc).*amcoeff).*Value_GSBF);
        
        D_set=[]; A_set=BS_index;  %A_set: active RRH set, D_set: inactive RRH set
        
        %% Process Deflation Procedure
        for A=0:L-2
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            %params.weight=(1/amcoeff)*ones(length(params.Active_index),1);
            params.weight=ones(length(params.Active_index),1); 
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if feasible==1
                
                TotalPower_GSBFHeuristic(lq)=(1/amcoeff)*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %recode current values
                TransmitPower_GSBFHeuristic(lq)=(1/amcoeff)*norm(Wsolution,'fro')^2;
                A_number_Heuristic(lq)=length(A_set);
                
                D_set=[D_set, BS_index(A+1)];  %%%process next selection: updata active RRH set and number
                A_set=setdiff([1:L],D_set); 
            else
                break;
            end
        end   
    end 
    
end


%% GSBF with Alternating Optimization Based Group Sparsity Design%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if GSBF_Alternating==true
    
    %% Compute Approximated Group Sparse Beamformer
    
    alternating_max=0; alternating_abs=10^99; 
    value1=10^99; value2=0;
    mu=1/L.*ones(L,1);
  
    while alternating_max<=10&alternating_abs>=10^(-3)
        
        alternating_max=alternating_max+1; %interation numbers
        value1=value2; %recode the old objective value
        
        params.Inactive_index=[]; 
        params.Active_index=[1:L];
        %params.weight=1./mu;
        params.weight=(1./mu).*(Pc(params.Active_index)./amcoeff(params.Active_index)); 
        params.rankone=false; %return general rank solution
        [Q_GSBF, feasible_GSBF] = powermin_cvx(params);
        
        if feasible_GSBF==false
            break;
        else
            value_temp=zeros(L,1);
            for l=1:L 
                temp1=0;
                for m=1:M
                    temp1=temp1+trace(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
                end
                value_temp(l)=temp1+10^(-3)/(M*sum(N_set));
            end
            mu=sqrt((Pc(params.Active_index)./amcoeff(params.Active_index)).*value_temp)./sum(sqrt((Pc(params.Active_index)./amcoeff(params.Active_index)).*value_temp)); %updata mu
            value2=((1./mu).*(Pc(params.Active_index)./amcoeff(params.Active_index)))'*value_temp; %new objective value
            alternating_abs=abs(value1-value2); %absolute value of the adjacent objective values
        end
    end
    
    if feasible_GSBF==false
        SU_channel_temp=0;  %This channel is infeasible
        break;
    else  
        %% RRH Ordering Based on the Approximated Group Sparse Beamformer
        Value_GSBF=zeros(L,1); 
        for l=1:L 
            temp1=0;
            for m=1:M
                temp1=temp1+trace(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
            end
            Value_GSBF(l)=temp1;
        end
        
        %% Sparsity Parttern Based RRH Ordering 
      %  [BS,BS_index]=sort(Value_GSBF);
        
        %% Proposed RRH Ordering
        ChannelGain=zeros(L,1);
        for l=1:L
            ChannelGain(l)=norm(vec(params.H(sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),:,:)),'fro')^2;
        end
        
        %[BS,BS_index]=sort(ChannelGain.*Value_GSBF);
        [BS,BS_index]=sort(((ChannelGain./Pc).*amcoeff).*Value_GSBF);
        
        %%
        D_set=[]; A_set=BS_index;  %A_set: active RRH set, D_set: inactive RRH set
        
        %% Process Deflation Procedure
        for A=0:L-2
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            %params.weight=(1/amcoeff)*ones(length(params.Active_index),1);
            params.weight=ones(length(params.Active_index),1); 
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if feasible==1
                
                TotalPower_GSBFAlternating(lq)=(1/amcoeff(1))*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %recode current values
                TransmitPower_GSBFAlternating(lq)=(1/amcoeff(1))*norm(Wsolution,'fro')^2;
                A_number_Alternating(lq)=length(A_set);
                
                D_set=[D_set, BS_index(A+1)];  %%%process next selection: updata active RRH set and number
                A_set=setdiff([1:L],D_set); 
                
            else
                break;
            end
        end   
    end 
    
end


%% Baseline Algorith with l1/l_infty Group Sparsity Design%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Baseline==true
    
    %% Compute Approximated Group Sparse Beamformer
    params.Inactive_index=[]; 
    params.Active_index=[1:L];
    %params.weight=(Pc(1)/amcoeff*length(params.Active_index))*ones(length(params.Active_index),1);
    params.weight=ones(length(params.Active_index),1);
    params.rankone=false; %return general rank solution
    [Q_GSBF, feasible_GSBF] = baseline_cvx(params);
    
    if feasible_GSBF==false
        SU_channel_temp=0;  %This channel is infeasible
        break;
    else  
        %% RRH Ordering Based on the Approximated Group Sparse Beamformer
        Value_GSBF=zeros(L,1); 
        for l=1:L 
            temp1=0;
            for m=1:M
                temp1=temp1+trace(Q_GSBF(sum(N_set(1:l-1))+1:sum(N_set(1:l)),sum(N_set(1:l-1))+1:sum(N_set(1:l)),m));
            end
            Value_GSBF(l)=temp1;
        end
        
        %% Sparsity Parttern Based RRH Ordering 
        [BS,BS_index]=sort(Value_GSBF);
        
        %% Proposed RRH Ordering
%         ChannelGain=zeros(L,1);
%         for l=1:L
%             ChannelGain(l)=norm(vec(params.H(sum(N_set(1:(l-1)))+1:sum(N_set(1:l)),:,:)),'fro')^2;
%         end
%         [BS,BS_index]=sort(ChannelGain.*Value_GSBF);
        
        %%
        
        D_set=[]; A_set=BS_index;  %A_set: active RRH set, D_set: inactive RRH set
        
        %% Process Deflation Procedure
        for A=0:L-2
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            %params.weight=(1/amcoeff)*ones(length(params.Active_index),1);
            params.weight=ones(length(params.Active_index),1);
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if feasible==1
                
                TotalPower_Baseline(lq)=(1/amcoeff(1))*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %recode current values
                TransmitPower_Baseline(lq)=(1/amcoeff(1))*norm(Wsolution,'fro')^2;
                A_number_Baseline(lq)=length(A_set);
                
                D_set=[D_set, BS_index(A+1)];  %%%process next selection: updata active RRH set and number
                A_set=setdiff([1:L],D_set); 
            else
                break;
            end
        end   
    end 
    
end

%% Exhaustive Search%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Exhaustive==true
    
    P_optimal=10^99; 
    D_set=[]; A_set=[1:L];  %A_set: active RRH set, D_set: inactive RRH set
    D_optimal=[]; A_optimal=[1:L];  % optimal active and inactive RRH set
    
    for A=0:L-1
        BS_select=nchoosek([1:L],A);
        feasible_temp=0; %recorder feasible times
        for s=1:nchoosek(L,A)
            
            A_set=setdiff([1:L],BS_select(s,:));   %A_set: active RRH set, D_set: inactive RRH set
            D_set=setdiff([1:L], A_set);
            
            params.Inactive_index=D_set; 
            params.Active_index=A_set;
            params.weight=(1/amcoeff)*ones(length(params.Active_index),1);
            params.rankone=true; %return rankone solution
            [Wsolution, feasible] = powermin_cvx(params); %power minimization given the active RRH set
            
            if  feasible==1   %%%if feasible
                temp1=(1/amcoeff)*norm(Wsolution,'fro')^2+sum(Pc(A_set));  %network power consumption
                feasible_temp=feasible_temp+1;
                if P_optimal>temp1
                    P_optimal=temp1;   %update optimal network power, active RRH set
                    D_optimal=D_set; A_optimal=A_set;  %record optimal set
                end
            end
        end
        
        if feasible_temp==0  %all the available sets are infeasible
            break;
        end
    end
    
%% Output for the Particular QoS
if P_optimal<10^50
    TotalPower_Exhaustive(lq)=P_optimal;  %Network Power consumption
    TransmitPower_Exhaustive(lq)=P_optimal-sum(Pc(A_optimal)); %Transmit Power Consumption
    A_number_Exhaustive(lq)=length(A_optimal);  %Number of active RRHs
else
    SU_channel_temp=0;  % this channel if infeasible
    break;
end
end


    end
    
%% Output for the Particular Channel Realization
if SU_channel_temp>0
    SU_counter=SU_counter+1;
    
    if CBF==true
        TotalPower_CBF_temp=TotalPower_CBF_temp+TotalPower_CBF;
        TransmitPower_CBF_temp=TransmitPower_CBF_temp+TransmitPower_CBF;
    end
    
    if GSBF_Heuristic==true
        TotalPower_GSBFHeuristic_temp=TotalPower_GSBFHeuristic_temp+TotalPower_GSBFHeuristic;
        TransmitPower_GSBFHeuristic_temp=TransmitPower_GSBFHeuristic_temp+TransmitPower_GSBFHeuristic;
        A_number_Heuristic_temp=A_number_Heuristic_temp+A_number_Heuristic;
    end
    
    if GSBF_Alternating==true
        TotalPower_GSBFAlternating_temp=TotalPower_GSBFAlternating_temp+TotalPower_GSBFAlternating;
        TransmitPower_GSBFAlternating_temp=TransmitPower_GSBFAlternating_temp+TransmitPower_GSBFAlternating; 
        A_number_Alternating_temp=A_number_Alternating_temp+A_number_Alternating;
    end
    
    if Baseline==true
        TotalPower_Baseline_temp=TotalPower_Baseline_temp+TotalPower_Baseline;
        TransmitPower_Baseline_temp=TransmitPower_Baseline_temp+TransmitPower_Baseline; 
        A_number_Baseline_temp=A_number_Baseline_temp+A_number_Baseline;
    end
    
    if Exhaustive==true
        TotalPower_Exhaustive_temp=TotalPower_Exhaustive_temp+TotalPower_Exhaustive;
        TransmitPower_Exhaustive_temp=TransmitPower_Exhaustive_temp+TransmitPower_Exhaustive; 
        A_number_Exhaustive_temp=A_number_Exhaustive_temp+A_number_Exhaustive;
    end
    
end
    
end

%% Final Results
if CBF==true
    TotalPower_CBF=TotalPower_CBF_temp./SU_counter;
    TransmitPower_CBF=TransmitPower_CBF_temp./SU_counter;
end

if GSBF_Heuristic==true
    TotalPower_GSBFHeuristic=TotalPower_GSBFHeuristic_temp./SU_counter;
    TransmitPower_GSBFHeuristic=TransmitPower_GSBFHeuristic_temp./SU_counter;
    A_number_Heuristic=A_number_Heuristic_temp/SU_counter;
end

if GSBF_Alternating==true
    TotalPower_GSBFAlternating=TotalPower_GSBFAlternating_temp./SU_counter;
    TransmitPower_GSBFAlternating=TransmitPower_GSBFAlternating_temp./SU_counter;
    A_number_Alternating=A_number_Alternating_temp/SU_counter;
end

if Baseline==true
    TotalPower_Baseline=TotalPower_Baseline_temp./SU_counter;
    TransmitPower_Baseline=TransmitPower_Baseline_temp./SU_counter;
    A_number_Baseline=A_number_Baseline_temp/SU_counter;
end

if Exhaustive==true
    TotalPower_Exhaustive=TotalPower_Exhaustive_temp./SU_counter;
    TransmitPower_Exhaustive=TransmitPower_Exhaustive_temp./SU_counter;
    A_number_Exhaustive=A_number_Exhaustive_temp/SU_counter;
end

%% Plot everage network power consumption%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
if CBF==true
    plot(Q,TotalPower_CBF,'b-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
     FronthaulPower_CBF=TotalPower_CBF-TransmitPower_CBF;
    hold on;
end

if GSBF_Heuristic==true
    plot(Q,TotalPower_GSBFHeuristic,'r-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    FronthaulPower_GSBFHeuristic=TotalPower_GSBFHeuristic-TransmitPower_GSBFHeuristic;
    hold on;
end

if GSBF_Alternating==true
    plot(Q,TotalPower_GSBFAlternating,'y-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    FronthaulPower_GSBFAlternating=TotalPower_GSBFAlternating-TransmitPower_GSBFAlternating;
    hold on;
end

if Baseline==true
    plot(Q,TotalPower_Baseline,'black-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    FronthaulPower_Baseline=TotalPower_Baseline-TransmitPower_Baseline;
    hold on;
end

if Exhaustive==true
    plot(Q,TotalPower_Exhaustive,'g-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    FronthaulPower_Exhaustive=TotalPower_Exhaustive-TransmitPower_Exhaustive;
    hold on;
end

h=legend('CBF', 'GSBF with Heuristic Optimization', 'GSBF with Alternating Optimization', 'Baseline', 'Exhaustive Search');  
xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Power Consumption [W]','fontsize',14,'fontweight','b','fontname','helvetica');


%% Plot everage transmit power consumption%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
if CBF==true
    plot(Q,TransmitPower_CBF,'b-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if GSBF_Heuristic==true
    plot(Q,TransmitPower_GSBFHeuristic,'r-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if GSBF_Alternating==true
    plot(Q,TransmitPower_GSBFAlternating,'y-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Baseline==true
    plot(Q,TransmitPower_Baseline,'black-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Exhaustive==true
    plot(Q,TransmitPower_Exhaustive,'g-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

h=legend('CBF', 'GSBF with Heuristic Optimization', 'GSBF with Alternating Optimization', 'Baseline', 'Exhaustive Search');  
xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Transmit Power Consumption [W]','fontsize',14,'fontweight','b','fontname','helvetica');


%% Plot everage number of active RRHs%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
if CBF==true
    plot(Q,L*ones(length(Q),1),'b-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if GSBF_Heuristic==true
    plot(Q,A_number_Heuristic,'r-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if GSBF_Alternating==true
    plot(Q,A_number_Alternating,'y-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Baseline==true
    plot(Q,A_number_Baseline,'black-o','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

if Exhaustive==true
    plot(Q,A_number_Exhaustive,'g-d','LineWidth',1.5, 'MarkerSize',10); %Network power consumptioin
    hold on;
end

h=legend('CBF', 'GSBF with Heuristic Optimization', 'GSBF with Alternating Optimization', 'Baseline', 'Exhaustive Search');  
xlabel('Target SINR [dB]','fontsize',12,'fontweight','b','fontname','helvetica');
ylabel('Average Number of Active RRHs','fontsize',14,'fontweight','b','fontname','helvetica');