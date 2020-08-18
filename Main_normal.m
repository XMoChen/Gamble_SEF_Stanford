clear all
Dir='E:\Data_Spike_LFP_Full\';
fileDir=dir(Dir);
fileDir_name={fileDir(3:end).name};
Dirfigure='D:\Projects\SEFCooling\figures_Normal\';

for i=[20:21,23:length(fileDir_name)];%[6:12,15:21,23:length(fileDir_name)] %13 14      22
    file_dir=[Dir,fileDir_name{i}];
    file=[file_dir,'\LFP.mat'];
    eventfile=[file_dir,'\Event.mat'];
    %      TimeFrequency_normal(file);
    %      TimeFrequencyPlot_normal(Dirfigure,fileDir_name{i},file,eventfile)
    
              Spect_normal(Dirfigure,fileDir_name{i},file,eventfile)
    %          SpectMI_normal(Dirfigure,fileDir_name{i},file,eventfile)
    %   Spect(Dirfigure,fileDir_name{i},file)
end


%%
TargDif=[];Targ_Sing=[];
ResultDif=[];ResultExpDif=[];

TargRaw=[];ResultRaw=[];ResultExpRaw=[];


TargInfo_SD=[];TargInfo_CD=[];
TargInfo_SD0=[];TargInfo_CD0=[];
TargInfo_SOpt=[];TargInfo_COpt=[];
TargInfo_SOpt0=[];TargInfo_COpt0=[];
ResultInfo_W=[];ResultInfo_R=[];
ResultInfo_W0=[];ResultInfo_R0=[];
ResultInfo_Rexp=[];ResultInfo_Rexp0=[];
ch_n=0;
TargMapMean=zeros(205,1024);
ResultMapMean=zeros(205,1024);
for i=[1:12,15:21,23:length(fileDir_name)]
    file_dir=[Dir,fileDir_name{i}];
    file=[file_dir,'\LFP.mat'];
    clear Spectrum SpectrumInfoResult SpectrumInfoTarg Map
    load(file,'Spectrum','SpectrumInfoResult','SpectrumInfoTarg','ChannelNumber','Map');
    %     Targ_Choice=cat(1,Targ_Choice,reshape(nanmean(Spectrum.ChoiceDirmean(:,1:4,:),2),[],size(Spectrum.ChoiceDirmean,3)));
    %     Targ_Sing=cat(1,Targ_Sing,reshape(nanmean(Spectrum.SingDirmean(:,1:4,:),2),[],size(Spectrum.ChoiceDirmean,3)));
    
    %     Result_all=cat(1,Result_all,reshape(nanmean(Spectrum.RewardMean(:,1:4,:),2),[],size(Spectrum.ChoiceDirmean,3))); %Spectrum.DifMean(:,3,:);%
    
    TargDif=cat(1,TargDif,reshape(Spectrum.DifMean(:,1,:),[],size(Spectrum.DifMean,3)));
    ResultDif=cat(1,ResultDif,reshape(Spectrum.DifMean(:,3,:),[],size(Spectrum.DifMean,3)));
    ResultExpDif=cat(1,ResultExpDif,reshape(Spectrum.DifMean(:,4,:),[],size(Spectrum.DifMean,3)));
    
    TargRaw=cat(1,TargRaw,reshape(Spectrum.RawMean(:,2,:),[],size(Spectrum.RawMean,3)));
    ResultRaw=cat(1,ResultRaw,reshape(Spectrum.RawMean(:,4,:),[],size(Spectrum.RawMean,3)));
    ResultExpRaw=cat(1,ResultExpRaw,reshape(Spectrum.RawMean(:,5,:),[],size(Spectrum.RawMean,3)));
    
    
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoTarg.Sig_SD);
    TargInfo_SD=cat(1,TargInfo_SD,SpectrumInfoTarg.Inf_RS_SD(:,:).*Sig_I);
    TargInfo_SD0=cat(1,TargInfo_SD0,SpectrumInfoTarg.Inf_RS_SD(:,:));
    
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoTarg.Sig_CD);
    TargInfo_CD=cat(1,TargInfo_CD,SpectrumInfoTarg.Inf_RS_CD(:,:).*Sig_I);
    TargInfo_CD0=cat(1,TargInfo_CD0,SpectrumInfoTarg.Inf_RS_CD(:,:));
    
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoTarg.Sig_SOpt);
    TargInfo_SOpt=cat(1,TargInfo_SOpt,SpectrumInfoTarg.Inf_RS_SOpt(:,:).*Sig_I);
    TargInfo_SOpt0=cat(1,TargInfo_SOpt0,SpectrumInfoTarg.Inf_RS_SOpt(:,:));
    
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoTarg.Sig_COpt);
    TargInfo_COpt=cat(1,TargInfo_COpt,SpectrumInfoTarg.Inf_RS_COpt(:,:).*Sig_I);
    TargInfo_COpt0=cat(1,TargInfo_COpt0,SpectrumInfoTarg.Inf_RS_COpt(:,:));
    
    
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoResult.Sig_W);
    ResultInfo_W=cat(1,ResultInfo_W,SpectrumInfoResult.Inf_RS_W(:,:).*Sig_I);
    ResultInfo_W0=cat(1,ResultInfo_W0,SpectrumInfoResult.Inf_RS_W(:,:));
    
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoResult.Sig_R);
    ResultInfo_R=cat(1,ResultInfo_R,SpectrumInfoResult.Inf_RS_R(:,:).*Sig_I);
    ResultInfo_R0=cat(1,ResultInfo_R0,SpectrumInfoResult.Inf_RS_R(:,:));
    
    
    clear Sig_I
    Sig_I=Sig_FDR(SpectrumInfoResult.Sig_Rexp);
    ResultInfo_Rexp=cat(1,ResultInfo_Rexp,SpectrumInfoResult.Inf_RS_Rexp(:,:).*Sig_I);
    ResultInfo_Rexp0=cat(1,ResultInfo_Rexp0,SpectrumInfoResult.Inf_RS_Rexp(:,:));
    
    
    for ch=1:ChannelNumber
        ch_n=ch_n+1;
        clear a b
        eval(['a=Map.TargMean',num2str(ch),';']);
        TargMapMean=TargMapMean+a;
        TargMapMean_d(i,:)=nanmean(a,2);
        
        eval(['b=Map.ResultMean',num2str(ch),';']);
        ResultMapMean=ResultMapMean+b;
        ResultMapMean_d(i,:)=nanmean(b,2);
        
    end
    
    
end



%%
f=linspace(0,500, size(TargRaw,2));

figure(1)
subplot(231)
plotstd(f,TargRaw,'b');hold on;
axis([5 150 -15 25]);
xlabel('Frequency (Hz)')
ylabel('Energy (dB)');
box off;set(gca,'TickDir','out');
title('Target')

subplot(232)
plotstd(f,ResultRaw,'b');hold on;
axis([5 150 -15 25]);
xlabel('Frequency (Hz)')
ylabel('Energy (dB)');
box off;set(gca,'TickDir','out');
title('Result')

subplot(233)
plotstd(f,ResultExpRaw,'b');hold on;
axis([5 150 -15 25]);
xlabel('Frequency (Hz)')
ylabel('Energy (dB)');
box off;set(gca,'TickDir','out');
title('Reward Expectation')

subplot(234)
plotstd(f,TargDif,'b');hold on;
axis([5 150 -1 1.5]);
xlabel('Frequency (Hz)')
ylabel('Energy (dB)');
box off;set(gca,'TickDir','out');

subplot(235)
plotstd(f,ResultDif,'b');hold on;
axis([5 150 -1 1.5]);
xlabel('Frequency (Hz)')
ylabel('Energy (dB)');
box off;set(gca,'TickDir','out');


subplot(236)
plotstd(f,ResultExpDif,'b');hold on;
axis([5 150 -1 1.5]);
xlabel('Frequency (Hz)')
ylabel('Energy (dB)');
box off;set(gca,'TickDir','out');
%%
f=linspace(1,500,size(Spectrum.TargetFFT,2));
f_d=downsample(f,2);
figure(2)
subplot(421)
plotstd(f_d(f_d<200),TargInfo_SD,'b');hold on;
title('Single Direction');box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)');
ylabel('Information (bit)');

subplot(422)
plotstd(f_d(f_d<200),TargInfo_CD,'b');hold on;
title('Choice Direction');box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)');
ylabel('Information (bit)');

subplot(423)
plotstd(f_d(f_d<200),TargInfo_SOpt,'b');hold on;
title('Single Option');box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)');
ylabel('Information (bit)');
axis([0 200 0 inf])

subplot(424)
plotstd(f_d(f_d<200),TargInfo_COpt,'b');hold on;
title('Choice Option');box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)');
ylabel('Information (bit)');

subplot(425)
plotstd(f_d(f_d<200),ResultInfo_W,'b');hold on;
title('Win/Loss');box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)');
ylabel('Information (bit)');

subplot(426)
plotstd(f_d(f_d<200),ResultInfo_W,'b');hold on;
title('Reward');box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)');
ylabel('Information (bit)');

subplot(427)
plotstd(f_d(f_d<200),ResultInfo_Rexp,'b');hold on;
title('Expectation');box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)');
ylabel('Information (bit)');

%%
figure(3)
subplot(221)
imagesc(f_d(f_d<200),1:size(TargInfo_SD,1),TargInfo_SD,[0 0.1]);
box off;set(gca,'TickDir','out');xlabel('Frequency (Hz)');
xlabel('Recordings');
colorbar;

subplot(222)
imagesc(f_d(f_d<200),1:size(TargInfo_SD,1),TargInfo_CD,[0 0.1]);
box off;set(gca,'TickDir','out');xlabel('Frequency (Hz)');
xlabel('Recordings');
colorbar;

subplot(223)
imagesc(f_d(f_d<200),1:size(TargInfo_SD,1),ResultInfo_Rexp,[0 0.2]);
box off;set(gca,'TickDir','out');xlabel('Frequency (Hz)');
xlabel('Recordings');
colorbar;

subplot(224)
imagesc(f_d(f_d<200),1:size(TargInfo_SD,1),ResultInfo_R,[0 0.2]);
box off;set(gca,'TickDir','out');xlabel('Frequency (Hz)');
xlabel('Recordings');
colorbar;

%%
f=f_d(f_d<200);
I_Theta=f>=4 & f<=8;
I_Alpha=f>=8 & f<=12;
I_Beta=f>=16 & f<=30;
I_Gamma=f>=30 & f<=80;
I_HGamma=f>=80 & f<=150;

Flag={'Alpha','Beta','Gamma','HGamma'};
figure(4)
for i=1:4
    subplot(2,2,i)
    eval(['I=I_',char(Flag(i)),';']);
    a=TargInfo_CD;
    a0=TargInfo_CD0;
    a1=hist(nanmean(a(:,I),2),0:0.01:0.15);
    a2=hist(nanmean(a0(:,I),2),0:0.01:0.15);
    
    bar(0:0.01:0.15,a2,'BarWidth',0.8, 'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
    bar(0:0.01:0.15,a1,'BarWidth',0.8, 'FaceColor',[0.2 0.2 0.2],'EdgeColor',[1 1 1]);
    
    axis([-0.02 0.15 0 100]);
    [~,p]=ttest(nanmean(a(:,I),2));
    title([char(Flag(i)),num2str(p,3)]);
    box off;set(gca,'TickDir','out');
    
end
%%
figure(5)
for i=1:4
    subplot(2,2,i)
    eval(['I=I_',char(Flag(i)),';']);
    a=ResultInfo_W;
    a0=ResultInfo_W0;
    a1=hist(nanmean(a(:,I),2),0:0.01:0.15);
    a2=hist(nanmean(a0(:,I),2),0:0.01:0.15);
    
    bar(0:0.01:0.15,a2,'BarWidth',0.8, 'FaceColor',[0.5 0.5 0.5],'EdgeColor',[1 1 1]);
    bar(0:0.01:0.15,a1,'BarWidth',0.8, 'FaceColor',[0.2 0.2 0.2],'EdgeColor',[1 1 1]);
    axis([-0.02 0.5 0 100]);
    [~,p]=ttest(nanmean(a(:,I),2));
    title([char(Flag(i)),num2str(p,3)]);
    box off;set(gca,'TickDir','out');
    
end


%%
figure(6)
t_MP=Map.t;
f_MP=Map.f;
clim=[-1 1];
subplot(2,2,1)
a=TargMapMean/ch_n;
tfplot(t_MP,f_MP,a,clim)
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-100 400 0 150]);
xlabel('Time from target onset (ms)');
ylabel('Frequency (Hz)');    colorbar;

% figure(gcf);set(gca,'YDir','norm');
subplot(222)
a=ResultMapMean/ch_n;
tfplot(t_MP,f_MP,a,clim)
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-300 400 0 150]);
xlabel('Time from result onset (ms)');
ylabel('Frequency (Hz)');    colorbar;

subplot(223)
h1=plot(f_MP,nanmean(TargMapMean_d,1),'b');hold on;
plotstd(f_MP,TargMapMean_d,'b');hold on;
axis([5 150 -1.5 1.5]);
box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');
title('Target') ;   

subplot(224)
h1=plot(f_MP,nanmean(ResultMapMean_d,1),'b');hold on;
plotstd(f_MP,ResultMapMean_d,'b');hold on;
axis([5 150 -1.5 1.5]);
box off;set(gca,'TickDir','out');
xlabel('Frequency (Hz)'); ylabel('Energy (dB)');
title('Target') ;  

%%

for im=1:6
    FigHandle = figure(im);
    print( FigHandle, '-djpeg', [Dirfigure,'SEF_Control_Sum',num2str(im)]);
    print( FigHandle, '-depsc', [Dirfigure,'SEF_Control_Sum',num2str(im)]);
    savefig( FigHandle, [Dirfigure,'SEF_Control_Sum',num2str(im)]);
    
end