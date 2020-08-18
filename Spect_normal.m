function Spect_normal(Dirfigure,name,file,eventfile)
clear Spectrum
close all
Fs=1000;T=1/Fs;
L=512;
NFFT= 2^nextpow2(L);
load(file,'LFPTargAll','LFPResultAll','LFPMoveAll','ChannelNumber');
load(eventfile,'TrialInfo');
ch_num=ChannelNumber;
for ch=1:ch_num
  clear Raw* Target Move Result
        eval(['Target=LFPTargAll.Channel',num2str(ch),';']);
        eval(['Move=LFPMoveAll.Channel',num2str(ch),';']);
        eval(['Result=LFPResultAll.Channel',num2str(ch),';']);

        RawBg=squeeze(Target(:,101:400))';
        RawTg=squeeze(Target(:,501:800))';
        RawMv=squeeze(Move(:,351:650))';
        RawRs=squeeze(Result(:,501:800))';
        RawRs_exp=squeeze(Result(:,201:500))';

        
  [Spectrum.BGFFT(ch,:,:),f1]=((pmtm(RawBg,4,NFFT,Fs,'eigen')));     
  [Spectrum.TargetFFT(ch,:,:),f1]=((pmtm(RawTg,4,NFFT,Fs,'eigen')));
  [Spectrum.MoveFFT(ch,:,:),f1]=((pmtm(RawMv,4,NFFT,Fs,'eigen')));
  [Spectrum.ResultFFT(ch,:,:),f1]=(pmtm(RawRs,4,NFFT,Fs,'eigen'));
  [Spectrum.ResultFFT_exp(ch,:,:),f1]=(pmtm(RawRs_exp,4,NFFT,Fs,'eigen'));

end 


Spectrum.BGFFT=10*log10(abs(Spectrum.BGFFT));
Spectrum.TargetFFT=10*log10(abs(Spectrum.TargetFFT));
Spectrum.MoveFFT=10*log10(abs(Spectrum.MoveFFT));
Spectrum.ResultFFT=10*log10(abs(Spectrum.ResultFFT));
Spectrum.ResultFFT_exp=10*log10(abs(Spectrum.ResultFFT_exp));


Spectrum.BGFFT((Spectrum.BGFFT)<-20)=-20;
Spectrum.TargetFFT((Spectrum.TargetFFT)<-20)=-20;
Spectrum.MoveFFT((Spectrum.MoveFFT)<-20)=-20;
Spectrum.ResultFFT((Spectrum.ResultFFT)<-20)=-20;
Spectrum.ResultFFT_exp((Spectrum.ResultFFT_exp)<-20)=-20;

% Spectrum.BGFFT(isinf(Spectrum.BGFFT))=nan;
% Spectrum.TargetFFT(isinf(Spectrum.TargetFFT))=nan;
% Spectrum.MoveFFT(isinf(Spectrum.MoveFFT))=nan;
% Spectrum.ResultFFT(isinf(Spectrum.ResultFFT))=nan;
% Spectrum.ResultFFT_exp(isinf(Spectrum.ResultFFT_exp))=nan;

Spectrum.TargetFFTDif=Spectrum.TargetFFT-Spectrum.BGFFT;
Spectrum.MoveFFTDif=Spectrum.MoveFFT-Spectrum.BGFFT;
Spectrum.ResultFFTDif=Spectrum.ResultFFT-Spectrum.BGFFT;
Spectrum.ResultFFTDif_exp=Spectrum.ResultFFT_exp-Spectrum.BGFFT;

% plot(f1,squeeze(nanmean(nanmean(Spectrum.BGFFT,1),3)),'k');hold on;
% plot(f1,squeeze(nanmean(nanmean(Spectrum.TargetFFT,1),3)),'b');hold on;
% plot(f1,squeeze(nanmean(nanmean(Spectrum.MoveFFT,1),3)),'r');hold on;
% plot(f1,squeeze(nanmean(nanmean(Spectrum.ResultFFT,1),3)),'m');hold on;
a=squeeze(nanmean(nanmean(Spectrum.ResultFFT,1),3));
% axis([5 200 min(a(f1<200))*1.2 max(a(f1<200))*1.2]);


%1-4 control 5-8 inactivation
for ch=1:ch_num
Spectrum.RawMean(ch,1,:)=squeeze(nanmean(nanmean(Spectrum.BGFFT(ch,:,:),1),3));
Spectrum.RawMean(ch,2,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFT(ch,:,:),1),3));
Spectrum.RawMean(ch,3,:)=squeeze(nanmean(nanmean(Spectrum.MoveFFT(ch,:,:),1),3));
Spectrum.RawMean(ch,4,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFT(ch,:,:),1),3));
Spectrum.RawMean(ch,5,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFT_exp(ch,:,:),1),3));

Spectrum.DifMean(ch,1,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFTDif(ch,:,:),1),3));
Spectrum.DifMean(ch,2,:)=squeeze(nanmean(nanmean(Spectrum.MoveFFTDif(ch,:,:),1),3));
Spectrum.DifMean(ch,3,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFTDif(ch,:,:),1),3));
Spectrum.DifMean(ch,4,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFTDif_exp(ch,:,:),1),3));



for dir=1:4
    clear I_dir
    I_dir=TrialInfo{:,15}==dir & TrialInfo{:,6}~=0;
    Spectrum.ChoiceDirmean(ch,dir,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFT(ch,:,I_dir),1),3));

    clear I_dir
    I_dir=TrialInfo{:,15}==dir & TrialInfo{:,6}==0;
    Spectrum.SingDirmean(ch,dir,:)=squeeze(nanmean(nanmean(Spectrum.TargetFFT(ch,:,I_dir),1),3));
   
end


%%
reward=[1,3,5,9];
for dir=1:4
    r=reward(dir);
    clear I_dir
    I_dir=TrialInfo{:,17}==r ;
    Spectrum.RewardMean(ch,dir,:)=squeeze(nanmean(nanmean(Spectrum.ResultFFT(ch,:,I_dir),1),3));

end

% f0=f1(f1<200)';
% figure(ch)
% subplot(221)
% plotstd(f0,squeeze(((Spectrum.TargetFFT(ch,(f1<200),:))))','b');hold on;
%  axis([5 200 min(a(f1<200))*1.2 max(a(f1<200))*1.6]);
% title('Target Ch1')
% subplot(222)
% plotstd(f0,squeeze(((Spectrum.ResultFFT(ch,(f1<200),:))))','b');hold on;
% axis([5 200 min(a(f1<200))*1.2 max(a(f1<200))*1.6]);
% title('Result');
% subplot(223)
% plotstd(f0,squeeze(((Spectrum.TargetFFTDif(ch,(f1<200),:))))','b');hold on;
% axis([5 200 -1 2]);
% 
% subplot(224)
% plotstd(f0,squeeze(((Spectrum.ResultFFTDif(ch,(f1<200),:))))','b');hold on;
% axis([5 200 -1 2]);

end


%%

figure(1)
for ch=1:ch_num

subplot(4,2,(ch-1)*2+1)
plot(f1,squeeze(Spectrum.ChoiceDirmean(ch,1:4,:)));
axis([5 200 -inf inf]);
title(['Choice ch',num2str(ch)]);

subplot(4,2,(ch-1)*2+2)
plot(f1,squeeze(Spectrum.SingDirmean(ch,1:4,:)));
axis([5 200 -inf inf]);
title('Single ');
end



figure(2)
for ch=1:ch_num
subplot(1,4,ch)
plot(f1,squeeze(Spectrum.RewardMean(ch,1:4,:)));
axis([5 200 -inf inf]);
title('Result');
end

d=strfind(name,'.')-1;

h=figure(1)
print( h, '-djpeg', [Dirfigure,name,'fft_direction']);
h=figure(2)
print( h, '-djpeg', [Dirfigure,name,'fft_result']);

% 
 save(file,'Spectrum','f1','-append');