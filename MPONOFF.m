function [rEnergyAll]=MPONOFF(X)
folderName ='data1/';
tag = 'test/';

Fs = 1000; % Sampling frequency
L  = 1024; % signal length

L=size(X,1);
t  = ((0:L-1))/Fs;
signalRange = [1 L]; % full range
importData(X,folderName,tag,signalRange,Fs);

% perform Gabor decomposition
Numb_points = L; % length of the signal
Max_iterations = 100; % number of iterations
runGabor(folderName,tag,Numb_points, Max_iterations);

f_cal = 0:Fs/L:Fs/2;
% f = 3:3:150;
% t=3:3:420;


for ch=1:size(X,3)
clear gaborInfo rEnergy_sample
gaborInfo = getGaborData(folderName,tag,ch);
for trial=1:size(X,2)  
trialNum=trial; % plot Trial number 
% Reconstruct signal
wrap=1;
atomList=[]; % all atoms

% if isempty(atomList)
%     disp(['Reconstructing trial ' num2str(trialNum) ', all atoms']);
% else
%     disp(['Reconstructing trial ' num2str(trialNum) ', atoms ' num2str(atomList(1)) ':' num2str(atomList(end))]);
% end

% reconstruct energy
clear rEnergy A B gaborData rEnergy_sub
gaborData=gaborInfo{trialNum}.gaborData;
rEnergy = reconstructEnergyFromAtomsMPP(gaborInfo{trialNum}.gaborData,L,wrap,atomList);
%rEnergy_base0=nanmean((rEnergy(:,20:(end-20))),2)*ones(1,512);
A=10*log10(abs(rEnergy));A(A<-20)=-20;
%B=nanmean(A(:,40:(80)),2)*ones(1,512);
% rEnergy_base(ch,trialNum,:,:)=B;
rEnergy_sub=A(f_cal<200,:);%-B;%-rEnergy_base;


%rEnergy_sub=rEnergy_sub(f_cal<200,20:(end-20));
% rEnergy_sample0=downsample(rEnergy_sub,2);
% rEnergy_sample(trialNum,:,:)=downsample(rEnergy_sample0',3)';
rEnergy_sample(trialNum,:,:)=rEnergy_sub;
end
rEnergyAll{ch}=rEnergy_sample;
end

