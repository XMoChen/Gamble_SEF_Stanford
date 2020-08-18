function TimeFrequency_normal(file)
%load(file,'Target','Move','Result','Infortable');
% X(1:2,:,:)=Target(:,:,280:791);
% X(3:4,:,:)=Result(:,:,280:791);
load(file,'LFPTargAll','LFPResultAll','LFPMoveAll','ChannelNumber');
%load(eventfile,'TrialInfo');
ch_num=ChannelNumber;
for ch=1:ch_num
  clear Raw* Target Move Result
        eval(['Target=LFPTargAll.Channel',num2str(ch),';']);
        eval(['Result=LFPResultAll.Channel',num2str(ch),';']);

        
        X(ch,:,:)=Target(:,:);
        X(ch+ch_num,:,:)=Result(:,:);
       
end 



Fs=1000;T=1/Fs;
L=1024;
f_cal = 0:Fs/L:Fs/2;
% f0=f_cal(1:90);
t_MP=(1:1024)-500;
f_MP=f_cal(f_cal<200);

X=permute(X,[ 3 2 1]);
[rEnergyAll]=MPCoolingProb(X);
% [rEnergyAll,f,t]=SpectrogramProb(X);
% t_MP=t;
% f_MP=f;

%%
for ch=1:ch_num
    clear A B
A=rEnergyAll{ch}; %A(A<-10)=-10;
B=repmat(squeeze(nanmean(A(:,:,t_MP>-300 & t_MP<-100),3)),[1 1 size(A,3)] );
Targ0=A-B;
% A=rEnergyAll{3};
% B=repmat(squeeze(nanmean(A(:,:,t_MP>-100 & t_MP<0),3)),[1 1 size(t_MP,2)] );
clear B0
B0=rEnergyAll{ch+ch_num};%B0(B0<-10)=-10;
% B=repmat(squeeze(nanmean(B0(:,:,t_MP>-300 & t_MP<-100),3)),[1 1 size(A,3)] );
Result0=B0-B;

eval(['Map.TargMean',num2str(ch),'=squeeze(nanmean(single(Targ0),1));'])
eval(['Map.ResultMean',num2str(ch),'=squeeze(nanmean(single(Result0),1));'])    


Targ0=Targ0(:,:,t_MP>-100 & t_MP<400);
Result0=Result0(:,:,t_MP>-400 & t_MP<400);
Targ_t=downsample(t_MP(t_MP>-100 & t_MP<400),10);
Result_t=downsample(t_MP(t_MP>-400 & t_MP<400),10);
Targ_f=downsample(f_MP,2);
Result_f=downsample(f_MP,2);

a=downsample3D(Targ0,1,2,10);
b=downsample3D(Result0,1,2,10);

eval(['Map.Targ',num2str(ch),'=single(a);'])
eval(['Map.Result',num2str(ch),'=single(b);'])

end
% clear A B
% A=rEnergyAll{2};%A(A<-2)=-2;
% B=repmat(squeeze(nanmean(A(:,:,t_MP>-100 & t_MP<0),3)),[1 1 size(t_MP,2)] );
% Map.Targ2=rEnergyAll{2}-B;
% B0=rEnergyAll{4};%B0(B0<-2)=-2;
% Map.Result2=B0-B;
%%
clim=[-1 1];
figure(1)
for ch=1:ch_num
subplot(4,2,1+(ch-1)*2)
eval(['a=Map.TargMean',num2str(ch),';']);
tfplot(t_MP,f_MP,a,clim)
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-100 400 0 200]);
% figure(gcf);set(gca,'YDir','norm');
subplot(4,2,2+(ch-1)*2)
eval(['a=Map.ResultMean',num2str(ch),';']);
tfplot(t_MP,f_MP,a,clim)
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-300 400 0 200]);
end



Map.t=t_MP;
Map.f=f_MP;
Map.Targ_t=Targ_t;
Map.Targ_f=Targ_f;
Map.Result_t=Result_t;
Map.Result_f=Result_f;
% 
save(file,'Map','-append');
% h=figure(1)
% print( h, '-djpeg', [Dirfigure,name,'fft_direction']);
end
 function Y=downsample3D(X,a1,a2,a3)
 
X=downsample(X,a1); 
X=permute(X,[2 1 3]);
X=downsample(X,a2);
X=permute(X,[3 1 2]);
X=downsample(X,a3);
Y=permute(X,[3 2 1]);
 end
% save(