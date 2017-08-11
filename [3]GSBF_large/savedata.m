Q=[8:-2:0]';  %QoS in dB

save('TotalPower_CBF.mat','TotalPower_CBF');
TotalPower_CBF_large=zeros(length(Q),2);
TotalPower_CBF_large(:,1)=Q;
TotalPower_CBF_large(:,2)=TotalPower_CBF;
save('TotalPower_CBF_large.dat','TotalPower_CBF_large','-ascii');

save('TotalPower_GSBFAlternating.mat','TotalPower_GSBFAlternating');
TotalPower_GSBFAlternating_large=zeros(length(Q),2);
TotalPower_GSBFAlternating_large(:,1)=Q;
TotalPower_GSBFAlternating_large(:,2)=TotalPower_GSBFAlternating;
save('TotalPower_GSBFAlternating_large.dat','TotalPower_GSBFAlternating_large','-ascii');

save('TotalPower_Baseline.mat','TotalPower_Baseline');
TotalPower_Baseline_large=zeros(length(Q),2);
TotalPower_Baseline_large(:,1)=Q;
TotalPower_Baseline_large(:,2)=TotalPower_Baseline;
save('TotalPower_Baseline_large.dat','TotalPower_Baseline_large','-ascii');


% save('TotalPower_Exhaustive.mat','TotalPower_Exhaustive');
% TotalPower_Exhaustive_data=zeros(length(Q),2);
% TotalPower_Exhaustive_data(:,1)=Q;
% TotalPower_Exhaustive_data(:,2)=TotalPower_Exhaustive;
% save('TotalPower_Exhaustive_data.dat','TotalPower_Exhaustive_data','-ascii');


%%
save('TransmitPower_CBF.mat','TransmitPower_CBF');
save('TransmitPower_GSBFAlternating.mat','TransmitPower_GSBFAlternating');
save('TransmitPower_Baseline.mat','TransmitPower_Baseline');
%save('TransmitPower_Exhaustive.mat','TransmitPower_Exhaustive');

save('FronthaulPower_CBF.mat','FronthaulPower_CBF');
save('FronthaulPower_GSBFAlternating.mat','FronthaulPower_GSBFAlternating');
save('FronthaulPower_Baseline.mat','FronthaulPower_Baseline');
%save('FronthaulPower_Exhaustive.mat','FronthaulPower_Exhaustive');

%%
save('A_number_Alternating.mat','A_number_Alternating');
save('A_number_Baseline.mat','A_number_Baseline');
%save('A_number_Exhaustive.mat','A_number_Exhaustive');