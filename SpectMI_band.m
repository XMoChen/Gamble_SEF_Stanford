function SpectMI_band(Dirfigure,name,file)
     clear Spectrum
     close all
load(file,'Spectrum','Infortable');

f=linspace(1,500,size(Spectrum.TargetFFT,2));
f_d=downsample(f,2);


A=Spectrum.TargetFFT; A=permute(A,[  2 1 3]);
A0=downsample(A,2);
A0=permute(A0,[2 1 3]);
A0=A0(:,f_d<200,:);

I_single=table2array(Infortable(:,8))==0;
I_choice=table2array(Infortable(:,8))~=0;
%%%% choice direction
s=table2array(Infortable(:,17));I_s=(s>0);
I_norm=table2array(Infortable(:,1))>35;
I_cool=table2array(Infortable(:,1))<10;


%I_Theta=f>=4 & f<=8;
I_Alpha=f>=2 & f<=12;
I_Beta=f>=16 & f<=30;
I_Gamma=f>=30 & f<=60;
I_HGamma=f>=80 & f<=150;

Flag={'Alpha','Beta','Gamma','HGamma'};

%%%%% For direction
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for b=1:4
       eval(['I=',Flag(b),';']);
       a=nanmean(A0_c(I,:),1);
       
       [Inf_RS_norm_SD_b(ch,b),Sig0_norm_SD(ch,b)]=MutualInfo_NM(a(I_s & I_norm & I_single),s(I_s & I_norm & I_single),4,4);
       [Inf_RS_inact_SD_b(ch,b),Sig0_inact_SD(ch,b)]=MutualInfo_NM(a(I_s & I_cool & I_single),s(I_s & I_cool & I_single),4,4);

       [Inf_RS_norm_CD_b(ch,b),Sig0_norm_CD(ch,b)]=MutualInfo_NM(a(I_s & I_norm & I_choice),s(I_s & I_norm & I_choice),4,4);
       [Inf_RS_inact_CD_b(ch,b),Sig0_inact_CD(ch,b)]=MutualInfo_NM(a(I_s & I_cool & I_choice),s(I_s & I_cool & I_choice),4,4);


   end
end
SpectrumInfoTarg.Inf_RS_norm_SD_band=Inf_RS_norm_SD_b;
SpectrumInfoTarg.Inf_RS_inact_SD_band=Inf_RS_inact_SD_b;
SpectrumInfoTarg.Inf_RS_norm_CD_band=Inf_RS_norm_CD_b;
SpectrumInfoTarg.Inf_RS_inact_CD_band=Inf_RS_inact_CD_b;

% SpectrumInfoTarg.Inf_RS_norm_Opt=Inf_RS_norm_Opt;
% SpectrumInfoTarg.Inf_RS_inact_Opt=Inf_RS_inact_Opt;
SpectrumInfoTarg.Sig_norm_SD_band=Sig0_norm_SD_b;
SpectrumInfoTarg.Sig_inact_SD_band=Sig0_inact_SD_b;
SpectrumInfoTarg.Sig_norm_CD_band=Sig0_norm_CD_b;
SpectrumInfoTarg.Sig_inact_CD_band=Sig0_inact_CD_b;

% SpectrumInfoTarg.Sig_norm_Opt=Sig0_norm_Opt;
% SpectrumInfoTarg.Sig_inact_Opt=Sig0_inact_Opt;
% %%%%% For Option
% s=table2array(Infortable(:,9));I_s=(s>0);
% ss=(s>1)+2*(s==1);
% for ch=1:size(Spectrum.TargetFFT)
%    A0_c=squeeze( A0(ch,:,:));
%    for ff=1:size(A0_c,1)
%        a=A0_c(ff,:);
%        [Inf_RS_norm_Opt(ch,ff),Sig0_norm_Opt(ch,ff)]=MutualInfo_NM(a(I_s & I_norm),ss(I_s & I_norm),7,7);
%        [Inf_RS_inact_Opt(ch,ff),Sig0_inact_Opt(ch,ff)]=MutualInfo_NM(a(I_s & I_cool),ss(I_s & I_cool),7,7);
% 
%    end
% end       
figure(1)
subplot(241)
bar([Inf_RS_norm_SD_b(1,:);Inf_RS_inact_SD_b(1,:)]);hold on;
title('ch1 SD')
axis([0 200 0 0.2]);

subplot(242)
bar([Inf_RS_norm_SD_b(1,:);Inf_RS_inact_SD_b(2,:)]);hold on;
axis([0 200 0 0.2]);
title('ch2 SD')

subplot(243)
bar([Inf_RS_norm_CD_b(1,:);Inf_RS_inact_CD_b(1,:)]);hold on;
title('ch1 CD')
axis([0 200 0 0.2]);

subplot(244)
bar([Inf_RS_norm_CD_b(2,:);Inf_RS_inact_CD_b(2,:)]);hold on;
axis([0 200 0 0.2]);
title('ch2 CD')


%%
%%%%%%%%%%%%%%%%%%%%%%%% Result %%%%%%%%%%%%%%%%%%%%%%%%%%
clear A A0 a
A=Spectrum.ResultFFT; A=permute(A,[  2 1 3]);
A0=downsample(A,2);
A0=permute(A0,[2 1 3]);
A0=A0(:,f_d<200,:);


clear Inf_RS*
%%%%  For result   Option/Option
s=table2array(Infortable(:,20));I_s=(s>0);
%ss=(s>1)+2*(s==1);
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for b=1:4
       eval(['I=',Flag(b),';']);
       a=nanmean(A0_c(I,:),1);
       [Inf_RS_norm_R_b(ch,b),Sig0_norm_R_b(ch,b)]=MutualInfo_NM(a(I_s & I_norm),s(I_s & I_norm),4,4);
       [Inf_RS_inact_R_b(ch,b),Sig0_inact_R_b(ch,b)]=MutualInfo_NM(a(I_s & I_cool),s(I_s & I_cool),4,4);

   end
end


%%%%  For result   Win/Loss
s=table2array(Infortable(:,20));I_s=(s>0);
ss=(s>1)+2*(s==1);
%%%%% For win/loss
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for b=1:4
       eval(['I=',Flag(b),';']);
       a=nanmean(A0_c(I,:),1);
       [Inf_RS_norm_W_b(ch,b),Sig0_norm_W_b(ch,b)]=MutualInfo_NM(a(I_s & I_norm),ss(I_s & I_norm),2,2);
       [Inf_RS_inact_W_b(ch,b),Sig0_inact_W_b(ch,b)]=MutualInfo_NM(a(I_s & I_cool),ss(I_s & I_cool),2,2);

   end
end

subplot(245)
bar([Inf_RS_norm_R_b(1,:);Inf_RS_inact_R_b(1,:)]);hold on;
legend({'N','In'});
title('ch1 Win/Loss')

subplot(246)
bar([Inf_RS_norm_R_b(2,:);Inf_RS_inact_R_b(2,:)]);hold on;
title('ch2 Win/Loss')

subplot(247)
bar([Inf_RS_norm_W_b(1,:);Inf_RS_inact_W_b(1,:)]);hold on;
title('ch1 Opt')

subplot(248)
bar([Inf_RS_norm_R_b(2,:);Inf_RS_inact_R_b(2,:)]);hold on;
title('ch2 Opts')


SpectrumInfoResult.Inf_RS_norm_W_band=Inf_RS_norm_W_b;
SpectrumInfoResult.Inf_RS_inact_W_band=Inf_RS_inact_W_b;
SpectrumInfoResult.Inf_RS_norm_R_band=Inf_RS_norm_R_b;
SpectrumInfoResult.Inf_RS_inact_R_band=Inf_RS_inact_R_b;
SpectrumInfoResult.Sig_norm_W_band=Sig0_norm_W_b;
SpectrumInfoResult.Sig_inact_W_band=Sig0_inact_W_b;
SpectrumInfoResult.Sig_norm_R_band=Sig0_norm_R_b;
SpectrumInfoResult.Sig_inact_R_band=Sig0_inact_R_b;

d=strfind(name,'.')-1;
h=figure(1)
print( h, '-djpeg', [Dirfigure,name(1:d),'MI']);
% h=figure(2)
% print( h, '-djpeg', [Dirfigure,name(1:d),'MPCh2']);
% SpectrumInfo.Inf_RS_norm_D=Inf_RS_norm_D;
% SpectrumInfo.Inf_RS_inact_D=Inf_RS_inact_D;
save(file,'SpectrumInfoResult','SpectrumInfoTarg','-append');