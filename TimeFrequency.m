function TimeFrequency(file,Dirfigure,name)
load(file,'Target','Move','Result','Infortable');
% X(1:2,:,:)=cat(3,zeros(size(Target,1),size(Target,2),124),Target);
% X(3:4,:,:)=cat(3,zeros(size(Target,1),size(Target,2),124),Result);
X(1:2,:,:)=Target;
X(3:4,:,:)=Result;
Fs=1000;T=1/Fs;
L=1024;
f_cal = 0:Fs/L:Fs/2;
% f0=f_cal(1:90);
t_MP=(1:1024)-500;
f_MP=f_cal(f_cal<200);

X=permute(X,[ 3 2 1]);
[rEnergyAll]=MPCoolingProbAdd1(X);
% [rEnergyAll,f,t]=SpectrogramProb(X);
% t_MP=t;
% f_MP=f;
Control=Infortable{:,1}>=30 & Infortable{:,4}==0 & Infortable{:,5}==0 & Infortable{:,27}>0;
Inact=Infortable{:,1}<=15 & Infortable{:,4}==0 & Infortable{:,5}==0 & Infortable{:,27}>0;
%%
ch_num=size(Target,1);
for ch=1:ch_num
    clear A B
A=rEnergyAll{ch}; %A(A<-10)=-10;
B=repmat(squeeze(nanmean(A(:,:,t_MP>-300 & t_MP<-50),3)),[1 1 size(A,3)] );
Targ0=A-B;
% A=rEnergyAll{3};
% B=repmat(squeeze(nanmean(A(:,:,t_MP>-100 & t_MP<0),3)),[1 1 size(t_MP,2)] );
clear B0
B0=rEnergyAll{ch+ch_num};%B0(B0<-10)=-10;
% B=repmat(squeeze(nanmean(B0(:,:,t_MP>-300 & t_MP<-100),3)),[1 1 size(A,3)] );
Result0=B0-B;

eval(['Map.TargMean',num2str(ch),'_C=squeeze(nanmean(single(Targ0(Control,:,:)),1));'])
eval(['Map.TargMean',num2str(ch),'_I=squeeze(nanmean(single(Targ0(Inact,:,:)),1));'])

eval(['Map.ResultMean',num2str(ch),'_C=squeeze(nanmean(single(Result0(Control,:,:)),1));'])    
eval(['Map.ResultMean',num2str(ch),'_I=squeeze(nanmean(single(Result0(Inact,:,:)),1));'])    


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
clim=[-0.5 0.5];
figure(1)
for ch=1:ch_num
subplot(2,4,1+(ch-1)*4)
eval(['a=Map.TargMean',num2str(ch),'_C;']);
tfplot(t_MP,f_MP,a,clim)
title('Targ control')
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-100 400 0 200]);
subplot(2,4,2+(ch-1)*4)
eval(['a=Map.TargMean',num2str(ch),'_I;']);
tfplot(t_MP,f_MP,a,clim)
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
title('Targ Inactivation')
axis([-100 400 0 200]);
% figure(gcf);set(gca,'YDir','norm');
subplot(2,4,3+(ch-1)*4)
eval(['a=Map.ResultMean',num2str(ch),'_C;']);
tfplot(t_MP,f_MP,a,clim)
title('Result Control')
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-300 400 0 200]);
subplot(2,4,4+(ch-1)*4)
eval(['a=Map.ResultMean',num2str(ch),'_I;']);
tfplot(t_MP,f_MP,a,clim)
title('Result Inactivation')
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
h=figure(1)
print( h, '-djpeg', [Dirfigure,name(1:end-4)]);
%end

% sav
% Fs=1000;T=1/Fs;
% L=1024;
% f_cal = 0:Fs/L:Fs/2;
% f0=f_cal(1:90);
% t_MP=(1:1024)-500;
% f_MP=f_cal(f_cal<200);
% % t_MP=downsample(-100:(512-120-20),3);
% % f_MP=downsample(f0,2);
% 
% X=permute(X,[ 3 2 1]);
% [rEnergyAll]=MPCoolingProb(X);
% 
% 
% 
% A=rEnergyAll{1}; %A(A<-10)=-10;
% B=repmat(squeeze(nanmean(A(:,:,t_MP>-100 & t_MP<0),3)),[1 1 size(t_MP,2)] );
% Map.Targ1=rEnergyAll{1}-B;
% % A=rEnergyAll{3};
% % B=repmat(squeeze(nanmean(A(:,:,t_MP>-100 & t_MP<0),3)),[1 1 size(t_MP,2)] );
% B0=rEnergyAll{3};%B0(B0<-10)=-10;
% Map.Result1=B0-B;
% 
% clear A B
% A=rEnergyAll{2};%A(A<-2)=-2;
% B=repmat(squeeze(nanmean(A(:,:,t_MP>-100 & t_MP<0),3)),[1 1 size(t_MP,2)] );
% Map.Targ2=rEnergyAll{2}-B;
% B0=rEnergyAll{4};%B0(B0<-2)=-2;
% Map.Result2=B0-B;
% 
% figure()
% subplot(221)
% imagesc(squeeze(nanmean(Map.Targ1,1)))
% subplot(222)
% imagesc(squeeze(nanmean(Map.Result1,1)))
% subplot(223)
% imagesc(squeeze(nanmean(Map.Targ2,1)))
% subplot(224)
% imagesc(squeeze(nanmean(Map.Result2,1)))
% 
% Map.t=t_MP;
% Map.f=f_MP;
% 
% save(file,'Map','-append');
% 
% % save(