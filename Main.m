clear all
Dir='D:\Projects\SEFCooling\data\';
fileDir=dir(Dir);
fileDir_name={fileDir(4:25).name};
Dirfigure='D:\Projects\SEFCooling\figures\';

for i=1:length(fileDir_name)%[1:7,21:length(fileDir_name)]  %  9
    file=[Dir,fileDir_name{i}];
    %    AlignEvent(file);
            Spect(Dirfigure,fileDir_name{i},file);
    %      SpectMI(Dirfigure,fileDir_name{i},file)
    
    
   % TimeFrequency(file,Dirfigure,fileDir_name{i});
    %    TimeFrequencyPlot(Dirfigure,fileDir_name{i},file);
end
%%
% %%
TargMapMean_C=zeros(205,1024);
ResultMapMean_C=zeros(205,1024);
TargMapMean_I=zeros(205,1024);
ResultMapMean_I=zeros(205,1024);
ch_n=0;
for i=[1:length(fileDir_name)]
    file=[Dir,fileDir_name{i}];
    ChannelNumber=2;
    clear Spectrum
    load(file,'Spectrum','SpectrumInfoResult','SpectrumInfoTarg','Map');
    Targ_Choice_con((2*i-1):(2*i),:,:)=Spectrum.ChoiceDirmean(:,1:4,:);
    Targ_Choice_inact((2*i-1):(2*i),:,:)=Spectrum.ChoiceDirmean(:,5:8,:);
    Targ_Sing_con((2*i-1):(2*i),:,:)=Spectrum.SingDirmean(:,1:4,:);
    Targ_Sing_inact((2*i-1):(2*i),:,:)=Spectrum.SingDirmean(:,5:8,:);
    
    Result_con((2*i-1):(2*i),:,:)=Spectrum.RewardMean(:,1:4,:); %Spectrum.DifMean(:,3,:);%
    Result_inact((2*i-1):(2*i),:,:)=Spectrum.RewardMean(:,5:8,:); %Spectrum.DifMean(:,6,:);
    
    
    TargRaw_con((2*i-1):(2*i),:,:)=Spectrum.RawMean(:,2,:);
    TargRaw_inact((2*i-1):(2*i),:,:)=Spectrum.RawMean(:,6,:);
    ResultRaw_con((2*i-1):(2*i),:,:)=Spectrum.RawMean(:,4,:);
    ResultRaw_inact((2*i-1):(2*i),:,:)=Spectrum.RawMean(:,9,:);
    ResultExpRaw_con((2*i-1):(2*i),:,:)=Spectrum.RawMean(:,5,:);
    ResultExpRaw_inact((2*i-1):(2*i),:,:)=Spectrum.RawMean(:,10,:);
    
    
    TargDif_con((2*i-1):(2*i),:,:)=Spectrum.DifMean(:,1,:);
    TargDif_inact((2*i-1):(2*i),:,:)=Spectrum.DifMean(:,5,:);
    ResultDif_con((2*i-1):(2*i),:,:)=Spectrum.DifMean(:,3,:);
    ResultDif_inact((2*i-1):(2*i),:,:)=Spectrum.DifMean(:,7,:);
    ResultExpDif_con((2*i-1):(2*i),:,:)=Spectrum.DifMean(:,4,:);
    ResultExpDif_inact((2*i-1):(2*i),:,:)=Spectrum.DifMean(:,8,:);
    
    
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoTarg.Sig_norm_SD);
    TargInfo_SD_con((2*i-1):(2*i),:)=SpectrumInfoTarg.Inf_RS_norm_SD(:,:).*Sig_I;
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoTarg.Sig_inact_SD);
    TargInfo_SD_inact((2*i-1):(2*i),:)=SpectrumInfoTarg.Inf_RS_inact_SD(:,:).*Sig_I;
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoTarg.Sig_norm_CD);
    TargInfo_CD_con((2*i-1):(2*i),:)=SpectrumInfoTarg.Inf_RS_norm_CD(:,:).*Sig_I;
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoTarg.Sig_inact_CD);
    TargInfo_CD_inact((2*i-1):(2*i),:)=SpectrumInfoTarg.Inf_RS_inact_CD(:,:).*Sig_I;
    
    
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoResult.Sig_norm_W);
    ResultInfo_W_con((2*i-1):(2*i),:)=SpectrumInfoResult.Inf_RS_norm_W(:,:).*Sig_I;
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoResult.Sig_inact_W);
    ResultInfo_W_inact((2*i-1):(2*i),:)=SpectrumInfoResult.Inf_RS_inact_W(:,:).*Sig_I;
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoResult.Sig_norm_R);
    ResultInfo_R_con((2*i-1):(2*i),:)=SpectrumInfoResult.Inf_RS_norm_R(:,:).*Sig_I;
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoResult.Sig_inact_R);
    ResultInfo_R_inact((2*i-1):(2*i),:)=SpectrumInfoResult.Inf_RS_inact_R(:,:).*Sig_I;
    
    
    for ch=1:ChannelNumber
        ch_n=ch_n+1;
        clear a b
        eval(['a=Map.TargMean',num2str(ch),'_C;']);
        TargMapMean_C=TargMapMean_C+a;
        TargMapMean_C_d(i,:)=nanmean(a,2);
        
        eval(['b=Map.ResultMean',num2str(ch),'_C;']);
        ResultMapMean_C=ResultMapMean_C+b;
        ResultMapMean_C_d(i,:)=nanmean(b,2);
        
        clear a b
        eval(['a=Map.TargMean',num2str(ch),'_I;']);
        TargMapMean_I=TargMapMean_I+a;
        TargMapMean_I_d(i,:)=nanmean(a,2);
        
        eval(['b=Map.ResultMean',num2str(ch),'_I;']);
        ResultMapMean_I=ResultMapMean_I+b;
        ResultMapMean_I_d(i,:)=nanmean(b,2);
        
    end
    
end
%
%%
Targ_con_mean=squeeze(nanmean(cat(2,Targ_Choice_con,Targ_Sing_con),2));
Targ_inact_mean=squeeze(nanmean(cat(2,Targ_Choice_inact,Targ_Sing_inact),2));
Result_con_mean=squeeze(nanmean(Result_con,2));
Result_inact_mean=squeeze(nanmean(Result_inact,2));


TargRaw_con_mean=squeeze(nanmean(TargRaw_con,2));
TargRaw_inact_mean=squeeze(nanmean(TargRaw_inact,2));
ResultRaw_con_mean=squeeze(nanmean(ResultRaw_con,2));
ResultRaw_inact_mean=squeeze(nanmean(ResultRaw_inact,2));
ResultExpRaw_con_mean=squeeze(nanmean(ResultExpRaw_con,2));
ResultExpRaw_inact_mean=squeeze(nanmean(ResultExpRaw_inact,2));

TargDif_con_mean=squeeze(nanmean(TargDif_con,2));
TargDif_inact_mean=squeeze(nanmean(TargDif_inact,2));
ResultDif_con_mean=squeeze(nanmean(ResultDif_con,2));
ResultDif_inact_mean=squeeze(nanmean(ResultDif_inact,2));
ResultExpDif_con_mean=squeeze(nanmean(ResultExpDif_con,2));
ResultExpDif_inact_mean=squeeze(nanmean(ResultExpDif_inact,2));
f=linspace(0,500, size(Targ_con_mean,2));
%%
figure(1)
subplot(231)
h1=plot(f,nanmean(TargRaw_con_mean,1),'b');hold on;
h2=plot(f,nanmean(TargRaw_inact_mean),'r');hold on;
plotstd(f,TargRaw_con_mean,'b');hold on;
plotstd(f,TargRaw_inact_mean,'r');hold on;
axis([5 200 -15 inf]);
legend([h1 h2],{'Contr','Inact'});
box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');
title('Target')

subplot(232)
plotstd(f,ResultRaw_con_mean,'b');hold on;
plotstd(f,ResultRaw_inact_mean,'r');hold on;
axis([5 200 -15 inf]);
box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');
title('Reward')


subplot(233)
plotstd(f,ResultExpRaw_con_mean,'b');hold on;
plotstd(f,ResultExpRaw_inact_mean,'r');hold on;
axis([5 200 -15 inf]);
box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');
title('Reward Expectation')

subplot(234)
plotstd(f,TargDif_con_mean,'b');hold on;
plotstd(f,TargDif_inact_mean,'r');hold on;
axis([5 200 -inf inf]);
box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');

subplot(235)
plotstd(f,ResultDif_con_mean,'b');hold on;
plotstd(f,ResultDif_inact_mean,'r');hold on;
axis([5 200 -inf inf]);box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');

subplot(236)
plotstd(f,ResultExpDif_con_mean,'b');hold on;
plotstd(f,ResultExpDif_inact_mean,'r');hold on;
axis([5 200 -inf inf]);box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');
%
%%
f=linspace(1,500,size(Spectrum.TargetFFT,2));
f_d=downsample(f,2);
figure(2)
subplot(221)
plotstd(f_d(f_d<200),TargInfo_SD_con,'b');hold on;
plotstd(f_d(f_d<200),TargInfo_SD_inact,'r');hold on;
box off;set(gca,'TickDir','out');
axis([0 200 0 0.03]);
subplot(222)
plotstd(f_d(f_d<200),TargInfo_CD_con,'b');hold on;
plotstd(f_d(f_d<200),TargInfo_CD_inact,'r');hold on;
box off;set(gca,'TickDir','out');
subplot(223)
plotstd(f_d(f_d<200),ResultInfo_W_con,'b');hold on;
plotstd(f_d(f_d<200),ResultInfo_W_inact,'r');hold on;
box off;set(gca,'TickDir','out');
subplot(224)
plotstd(f_d(f_d<200),ResultInfo_R_con,'b');hold on;
plotstd(f_d(f_d<200),ResultInfo_R_inact,'r');hold on;
box off;set(gca,'TickDir','out');

I_Theta=f>=4 & f<=8;
I_Alpha=f>=2 & f<=12;
I_Beta=f>=16 & f<=30;
I_Gamma=f>=30 & f<=60;
I_HGamma=f>=80 & f<=150;

Flag={'Alpha','Beta','Gamma','HGamma'};
figure(3)
for i=1:4
    subplot(2,2,i)
    eval(['I=I_',char(Flag(i)),';']);
    a=Targ_con_mean-Targ_inact_mean;
    a0=hist(nanmean(a(:,I),2),-4:0.5:4);
    bar(-4:0.5:4,a0,'BarWidth',0.8, 'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
    axis([-4 4 0 35]);
    [~,p]=ttest(nanmean(a(:,I),2));
    title([char(Flag(i)),num2str(p,3)]);
    box off;set(gca,'TickDir','out');
    
end

figure(4)
for i=1:4
    subplot(2,2,i)
    eval(['I=I_',char(Flag(i)),';']);
    a=Result_con_mean-Result_inact_mean;
    a0=hist(nanmean(a(:,I),2),-4:0.5:4);
    bar(-4:0.5:4,a0,'BarWidth',0.8, 'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
    axis([-4 4 0 35]);
    [~,p]=ttest(nanmean(a(:,I),2));
    title([char(Flag(i)),num2str(p,3)]);
    box off;set(gca,'TickDir','out');
    
end


%%
figure(5)
t_MP=Map.t;
f_MP=Map.f;
clim=[-0.5 0.5];
subplot(2,2,1)
a=TargMapMean_C/ch_n;
tfplot(t_MP,f_MP,a,clim)
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-100 400 0 150]);
xlabel('Time from target onset (ms)');
ylabel('Frequency (Hz)');
title('Target Control');
colorbar;

% figure(gcf);set(gca,'YDir','norm');
subplot(222)
a=ResultMapMean_C/ch_n;
tfplot(t_MP,f_MP,a,[-0.4 0.4])
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-300 400 0 150]);
xlabel('Time from result onset (ms)');
ylabel('Frequency (Hz)');
title('Result Control');
colorbar;

%clim=[-0.1 0.1];

subplot(2,2,3)
a=TargMapMean_I/ch_n;
tfplot(t_MP,f_MP,a,clim)
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-100 400 0 150]);
xlabel('Time from target onset (ms)');
ylabel('Frequency (Hz)');    colorbar;
title('Result inactivation');

% figure(gcf);set(gca,'YDir','norm');
subplot(224)
a=ResultMapMean_I/ch_n;
tfplot(t_MP,f_MP,a,[min(a(:)) max(a(:))*0.1]);
colorbar;
%  imagesc(squeeze(nanmean(Map.Targ1,+1)))
axis([-300 400 0 150]);
xlabel('Time from result onset (ms)');
ylabel('Frequency (Hz)');
title('Result Inactivation');

%%
figure(6)
subplot(221)
h1=plot(f_MP,nanmean(TargMapMean_C_d,1),'b');hold on;
h2=plot(f_MP,nanmean(TargMapMean_I_d),'r');hold on;
plotstd(f_MP,TargMapMean_C_d,'b');hold on;
plotstd(f_MP,TargMapMean_I_d,'r');hold on;
axis([5 150 -5 1.5]);
legend([h1 h2],{'Contr','Inact'});
box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');
title('Target')
colorbar;

subplot(222)
h1=plot(f_MP,nanmean(ResultMapMean_C_d,1),'b');hold on;
h2=plot(f_MP,nanmean(ResultMapMean_I_d),'r');hold on;
plotstd(f_MP,ResultMapMean_C_d,'b');hold on;
plotstd(f_MP,ResultMapMean_I_d,'r');hold on;
axis([5 150 -5 1.5]);
legend([h1 h2],{'Contr','Inact'});
box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');
title('Result') ;   colorbar;

%%
% %%
% close all
% for c=1:size(Targ_Choice_con,1)
%     figure(c)
%     subplot(221)
%     plot(f,squeeze(Targ_Choice_con(c,:,:)));
%     a=Targ_Choice_con(c,:,:);
%     axis([0 200 min(a(:)) max(a(:))]);
%     subplot(222)
%     plot(f,squeeze(Targ_Choice_inact(c,:,:)));
%     axis([0 200 min(a(:)) max(a(:))]);
%     subplot(223)
%     a=Result_con(c,:,:);
%     plot(f,squeeze(Result_con(c,:,:)));
%     axis([0 200 min(a(:)) max(a(:))]);
%     subplot(224)
%     plot(f,squeeze(Result_inact(c,:,:)));
%     axis([0 200 min(a(:)) max(a(:))]);
% end
%
for im=1:6
    FigHandle = figure(im);
    print( FigHandle, '-djpeg', [Dirfigure,'SEF_CI_Sum',num2str(im)]);
    print( FigHandle, '-depsc', [Dirfigure,'SEF_CI_Sum',num2str(im)]);
    savefig( FigHandle, [Dirfigure,'SEF_CI_Sum',num2str(im)]);
    
end