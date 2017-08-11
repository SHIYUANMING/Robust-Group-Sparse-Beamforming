Q=[8:-2:0]';  %QoS in dB

save('TotalPower_CBF.mat','TotalPower_CBF');
TotalPower_CBF_data=zeros(length(Q),2);
TotalPower_CBF_data(:,1)=Q;
TotalPower_CBF_data(:,2)=TotalPower_CBF;
save('TotalPower_CBF_data.dat','TotalPower_CBF_data','-ascii');

save('TotalPower_GSBFAlternating.mat','TotalPower_GSBFAlternating');
TotalPower_GSBFAlternating_data=zeros(length(Q),2);
TotalPower_GSBFAlternating_data(:,1)=Q;
TotalPower_GSBFAlternating_data(:,2)=TotalPower_GSBFAlternating;
save('TotalPower_GSBFAlternating_data.dat','TotalPower_GSBFAlternating_data','-ascii');

save('TotalPower_Baseline.mat','TotalPower_Baseline');
TotalPower_Baseline_data=zeros(length(Q),2);
TotalPower_Baseline_data(:,1)=Q;
TotalPower_Baseline_data(:,2)=TotalPower_Baseline;
save('TotalPower_Baseline_data.dat','TotalPower_Baseline_data','-ascii');


save('TotalPower_Exhaustive.mat','TotalPower_Exhaustive');
TotalPower_Exhaustive_data=zeros(length(Q),2);
TotalPower_Exhaustive_data(:,1)=Q;
TotalPower_Exhaustive_data(:,2)=TotalPower_Exhaustive;
save('TotalPower_Exhaustive_data.dat','TotalPower_Exhaustive_data','-ascii');


%%
save('TransmitPower_CBF.mat','TransmitPower_CBF');
save('TransmitPower_GSBFAlternating.mat','TransmitPower_GSBFAlternating');
save('TransmitPower_Baseline.mat','TransmitPower_Baseline');
save('TransmitPower_Exhaustive.mat','TransmitPower_Exhaustive');

%%
save('A_number_Alternating.mat','A_number_Alternating');
save('A_number_Baseline.mat','A_number_Baseline');
save('A_number_Exhaustive.mat','A_number_Exhaustive');