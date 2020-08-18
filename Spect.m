function Spect(Dirfigure,name,file)
clear Spectrum
close all
Fs=1000;T=1/Fs;
L=512;
NFFT= 2^nextpow2(L);
load(file,'Target','Move','Result','Infortable');
  
for ch=1:2
  clear Raw*
        RawBg=squeeze(Target(ch,:,201:500))';
        RawTg=squeeze(Target(ch,:,501:800))';
        RawMv=squeeze(Move(ch,:,351:650))';
        RawRs=squeeze(Result(ch,:,501:800))';
        RawRs_exp=squeeze(Result(ch,:,201:500))';
        
  [Spectrum.BGFFT(ch,:,:),f1]=((pmtm(RawBg,4,NFFT,Fs,'eigen')));     
  [Spectrum.TargetFFT(ch,:,:),f1]=((pmtm(RawTg,4,NFFT,Fs,'eigen')));
  [Spectrum.MoveFFT(ch,:,:),f1]=((pmtm(RawMv,4,NFFT,Fs,'eigen')));
  [Spectrum.ResultFFT(ch,:,:),f1]=(pmtm(RawRs,4,NFFT,Fs,'eigen'));
  [Spectrum.ResultExpFFT(ch,:,:),f1]=(pmtm(RawRs_exp,4,NFFT,Fs,'eigen'));

end 

% A=10*log10(abs(Spectrum.BGFFT)+10^-20);A(A<-50)=-50;
Spectrum.BGFFT=EngergyBit(Spectrum.BGFFT);
Spectrum.TargetFFT=EngergyBit(Spectrum.TargetFFT);
Spectrum.MoveFFT=EngergyBit(Spectrum.MoveFFT);
Spectrum.ResultFFT=EngergyBit(Spectrum.ResultFFT);
Spectrum.ResultExpFFT=EngergyBit(Spectrum.ResultExpFFT);

% Spectrum.BGFFT(isinf(Spectrum.BGFFT))=nan;
% Spectrum.TargetFFT(isinf(Spectrum.TargetFFT))=nan;
% Spectrum.MoveFFT(isinf(Spectrum.MoveFFT))=nan;
% Spectrum.ResultFFT(isinf(Spectrum.ResultFFT))=nan;
% Spectrum.ResultExpFFT(isinf(Spectrum.ResultExpFFT))=nan;

Spectrum.TargetFFTDif=Spectrum.TargetFFT-Spectrum.BGFFT;
Spectrum.MoveFFTDif=Spectrum.MoveFFT-Spectrum.BGFFT;
Spectrum.ResultFFTDif=Spectrum.ResultFFT-Spectrum.BGFFT;
Spectrum.ResultExpFFTDif=Spectrum.ResultExpFFT-Spectrum.BGFFT;

% plot(f1,squeeze(nanmean(nanmean(Spectrum.BGFFT,1),3)),'k');hold on;
% plot(f1,squeeze(nanmean(nanmean(Spectrum.TargetFFT,1),3)),'b');hold on;
% plot(f1,squeeze(nanmean(nanmean(Spectrum.MoveFFT,1),3)),'r');hold on;
% plot(f1,squeeze(nanmean(nanmean(Spectrum.ResultFFT,1),3)),'m');hold on;

a=squeeze(nanmean(nanmean(Spectrum.ResultFFT,1),3));

% axis([5 200 min(a(f1<200))*1.2 max(a(f1<200))*1.2]);
%%
Control=Infortable{:,1}>=30 & Infortable{:,4}==0 & Infortable{:,5}==0 & Infortable{:,27}>0;
Inact=Infortable{:,1}<=15 & Infortable{:,4}==0 & Infortable{:,5}==0 & Infortable{:,27}>0;


%1-4 control 5-8 inactivation
for ch=1:2
Spectrum.RawMean(ch,1,:)=squeeze(nanmean(nanmean(Spectrum.BGFFT(ch,:,Control),1),3));
Spectrum.RawMean(ch,2,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFT(ch,:,Control),1),3));
Spectrum.RawMean(ch,3,:)=squeeze(nanmean(nanmean(Spectrum.MoveFFT(ch,:,Control),1),3));
Spectrum.RawMean(ch,4,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFT(ch,:,Control),1),3));
Spectrum.RawMean(ch,5,:)=squeeze(nanmean(nanmean(Spectrum.ResultExpFFT(ch,:,Control),1),3));

Spectrum.RawMean(ch,6,:)=squeeze(nanmean(nanmean(Spectrum.BGFFT(ch,:,Inact),1),3));
Spectrum.RawMean(ch,7,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFT(ch,:,Inact),1),3));
Spectrum.RawMean(ch,8,:)=squeeze(nanmean(nanmean(Spectrum.MoveFFT(ch,:,Inact),1),3));
Spectrum.RawMean(ch,9,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFT(ch,:,Inact),1),3));
Spectrum.RawMean(ch,10,:)=squeeze(nanmean(nanmean(Spectrum.ResultExpFFT(ch,:,:),1),3));



Spectrum.DifMean(ch,1,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(ch,:,Control),1),3));
Spectrum.DifMean(ch,2,:)=squeeze(nanmean(nanmean(Spectrum.MoveFFTDif(ch,:,Control),1),3));
Spectrum.DifMean(ch,3,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFTDif(ch,:,Control),1),3));
Spectrum.DifMean(ch,4,:)=squeeze(nanmean(nanmean(Spectrum.ResultExpFFTDif(ch,:,Control),1),3));

Spectrum.DifMean(ch,5,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(ch,:,Inact),1),3));
Spectrum.DifMean(ch,6,:)=squeeze(nanmean(nanmean(Spectrum.MoveFFTDif(ch,:,Inact),1),3));
Spectrum.DifMean(ch,7,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFTDif(ch,:,Inact),1),3));
Spectrum.DifMean(ch,8,:)=squeeze(nanmean(nanmean(Spectrum.ResultExpFFTDif(ch,:,Inact),1),3));
end

for dir=1:4
    clear I_dir
    I_dir=Control & Infortable{:,17}==dir & Infortable{:,8}~=0;
    Spectrum.ChoiceDirmean(1,dir,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(1,:,I_dir),1),3));
    Spectrum.ChoiceDirmean(2,dir,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(2,:,I_dir),1),3));
    clear I_dir
    I_dir=Inact & Infortable{:,17}==dir & Infortable{:,8}~=0;
    Spectrum.ChoiceDirmean(1,dir+4,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(1,:,I_dir),1),3));
    Spectrum.ChoiceDirmean(2,dir+4,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(2,:,I_dir),1),3));

    I_dir=Control & Infortable{:,17}==dir & Infortable{:,8}==0;
    Spectrum.SingDirmean(1,dir,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(1,:,I_dir),1),3));
    Spectrum.SingDirmean(2,dir,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(2,:,I_dir),1),3));
    clear I_dir
    I_dir=Inact & Infortable{:,17}==dir & Infortable{:,8}==0;
    Spectrum.SingDirmean(1,dir+4,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(1,:,I_dir),1),3));
    Spectrum.SingDirmean(2,dir+4,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(2,:,I_dir),1),3));

end
%%

f0=f1(f1<200)';
figure(1)
subplot(231)
plotstd(f0,squeeze(((Spectrum.TargetFFT(1,(f1<200),Inact))))','r');hold on;
plotstd(f0,squeeze(((Spectrum.TargetFFT(1,(f1<200),Control))))','b');hold on;
[5 200 min(a(f1<200))*1.2 max(a(f1<200))*1.6]
axis([5 200 -60 10 ]);
title('Target Ch1')
subplot(232)
plotstd(f0,squeeze(((Spectrum.ResultFFT(1,(f1<200),Control))))','b');hold on;
plotstd(f0,squeeze(((Spectrum.ResultFFT(1,(f1<200),Inact))))','r');hold on;
axis([5 200 -60 10 ]);
title('Result');
subplot(233)
plotstd(f0,squeeze(((Spectrum.ResultExpFFT(1,(f1<200),Control))))','b');hold on;
plotstd(f0,squeeze(((Spectrum.ResultExpFFT(1,(f1<200),Inact))))','r');hold on;
axis([5 200 -60 10 ]);
title('Result Exp');

subplot(234)
plotstd(f0,squeeze(((Spectrum.TargetFFTDif(1,(f1<200),Control))))','b');hold on;
plotstd(f0,squeeze(((Spectrum.TargetFFTDif(1,(f1<200),Inact))))','r');hold on;
axis([5 200 -1 2]);
subplot(235)
plotstd(f0,squeeze(((Spectrum.ResultFFTDif(1,(f1<200),Control))))','b');hold on;
plotstd(f0,squeeze(((Spectrum.ResultFFTDif(1,(f1<200),Inact))))','r');hold on;
axis([5 200 -1 2]);
subplot(236)
plotstd(f0,squeeze(((Spectrum.ResultExpFFTDif(1,(f1<200),Control))))','b');hold on;
plotstd(f0,squeeze(((Spectrum.ResultExpFFTDif(1,(f1<200),Inact))))','r');hold on;
axis([5 200 -1 2]);


figure(2)
subplot(221)
plotstd(f0,squeeze(((Spectrum.TargetFFT(2,(f1<200),Control))))','b');hold on;
plotstd(f0,squeeze(((Spectrum.TargetFFT(2,(f1<200),Inact))))','r');hold on;
axis([5 200 -60 10 ]);
title('Target Ch2');
subplot(222)
plotstd(f0,squeeze(((Spectrum.ResultFFT(2,(f1<200),Control))))','b');hold on;
plotstd(f0,squeeze(((Spectrum.ResultFFT(2,(f1<200),Inact))))','r');hold on;
axis([5 200 -60 10 ]);
title('Result');
subplot(223)
plotstd(f0,squeeze(((Spectrum.TargetFFTDif(2,(f1<200),Control))))','b');hold on;
plotstd(f0,squeeze(((Spectrum.TargetFFTDif(2,(f1<200),Inact))))','r');hold on;
axis([5 200 -1 2]);

subplot(224)
plotstd(f0,squeeze(((Spectrum.ResultFFTDif(2,(f1<200),Control))))','b');hold on;
plotstd(f0,squeeze(((Spectrum.ResultFFTDif(2,(f1<200),Inact))))','r');hold on;
axis([5 200 -1 2]);


figure(3)
subplot(221)
plot(f1,squeeze(Spectrum.ChoiceDirmean(1,1:4,:)));
axis([5 200 -2 6]);
title('Choice Control ch1');
subplot(223)
plot(f1,squeeze(Spectrum.ChoiceDirmean(1,5:8,:)));
axis([5 200 -2 6]);
title('Choice Inact');

subplot(222)
plot(f1,squeeze(Spectrum.SingDirmean(1,1:4,:)));
axis([5 200 -2 6]);
title('Single Control');

subplot(224)
plot(f1,squeeze(Spectrum.SingDirmean(1,5:8,:)));
axis([5 200 -2 6]);
title('Single Inact');

figure(4)
subplot(221)
plot(f1,squeeze(Spectrum.ChoiceDirmean(2,1:4,:)));
axis([5 200 -2 6]);
title('Choice Control ch2');
subplot(223)
plot(f1,squeeze(Spectrum.ChoiceDirmean(2,5:8,:)));
axis([5 200 -2 6]);
title('Choice Inact');

subplot(222)
plot(f1,squeeze(Spectrum.SingDirmean(2,1:4,:)));
axis([5 200 -2 6]);
title('Single Control');

subplot(224)
plot(f1,squeeze(Spectrum.SingDirmean(2,5:8,:)));
axis([5 200 -2 6]);
title('Single Inact');
%%
reward=[1,3,5,9];
for dir=1:4
    r=reward(dir);
    clear I_dir
    I_dir=Control & Infortable{:,20}==r & Infortable{:,19}~=0 ;
    Spectrum.RewardMean(1,dir,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFTDif(1,:,I_dir),1),3));
    Spectrum.RewardMean(2,dir,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFTDif(2,:,I_dir),1),3));
    clear I_dir
    I_dir=Inact & Infortable{:,20}==r & Infortable{:,19}~=0 ;
    Spectrum.RewardMean(1,dir+4,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFTDif(1,:,I_dir),1),3));
    Spectrum.RewardMean(2,dir+4,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFTDif(2,:,I_dir),1),3));

end

figure(5)
subplot(221)
plot(f1,squeeze(Spectrum.RewardMean(1,1:4,:)));
co=get(gca,'ColorOrder');
set(gca,'ColorOrder',[0 0 0;0.5 0 0;0.7 0 0; 1 0 0 ],'NextPlot','replacechildren');
plot(f1,squeeze(Spectrum.RewardMean(1,1:4,:)));
axis([5 200 -2 6]);
title('Choice Control ch1');
subplot(223)
plot(f1,squeeze(Spectrum.RewardMean(1,5:8,:)));
co=get(gca,'ColorOrder');
set(gca,'ColorOrder',[0 0 0;0.5 0 0;0.7 0 0; 1 0 0 ],'NextPlot','replacechildren');
plot(f1,squeeze(Spectrum.RewardMean(1,5:8,:)));
axis([5 200 -2 6]);
title('Choice Inact');
subplot(222)
plot(f1,squeeze(Spectrum.RewardMean(2,1:4,:)));
co=get(gca,'ColorOrder');
set(gca,'ColorOrder',[0 0 0;0.5 0 0;0.7 0 0; 1 0 0 ],'NextPlot','replacechildren');
plot(f1,squeeze(Spectrum.RewardMean(2,1:4,:)));
axis([5 200 -2 6]);
title('Choice Control ch2');
subplot(224)
plot(f1,squeeze(Spectrum.RewardMean(2,5:8,:)));
co=get(gca,'ColorOrder');
set(gca,'ColorOrder',[0 0 0;0.5 0 0;0.7 0 0; 1 0 0 ],'NextPlot','replacechildren');
plot(f1,squeeze(Spectrum.RewardMean(2,5:8,:)));
axis([5 200 -2 6]);
title('Choice Inact');

d=strfind(name,'.')-1;
h=figure(1)
print( h, '-djpeg', [Dirfigure,name(1:d),'fft_ch1']);
h=figure(2)
print( h, '-djpeg', [Dirfigure,name(1:d),'fft_ch2']);
h=figure(3)
print( h, '-djpeg', [Dirfigure,name(1:d),'Dir_ch1']);
h=figure(4)
print( h, '-djpeg', [Dirfigure,name(1:d),'Dir_ch2']);
h=figure(5)
print( h, '-djpeg', [Dirfigure,name(1:d),'result']);
% 
 save(file,'Spectrum','f1','-append');