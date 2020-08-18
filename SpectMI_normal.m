function SpectMI_normal(Dirfigure,name,file,eventfile)
     clear Spectrum
     close all
load(file,'Spectrum');
load(eventfile);

f=linspace(1,500,size(Spectrum.TargetFFT,2));
f_d=downsample(f,2);


A=Spectrum.TargetFFT; A=permute(A,[  2 1 3]);
A0=downsample(A,2);
A0=permute(A0,[2 1 3]);
A0=A0(:,f_d<200,:);

I_single=table2array(TrialInfo(:,6))==0;
I_choice=table2array(TrialInfo(:,6))~=0;

%%%% choice direction
s=table2array(TrialInfo(:,15));I_s=(s>0);


%%%%% For direction
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_SD(ch,ff),Sig0_SD(ch,ff)]=MutualInfo_NMS(a(I_s & I_single ),s(I_s & I_single),4,4);
       [Inf_RS_CD(ch,ff),Sig0_CD(ch,ff)]=MutualInfo_NMS(a(I_s & I_choice ),s(I_s & I_choice),4,4);

   end
end



%%%%% For Option
s=table2array(TrialInfo(:,14));I_s=(s>0);
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_SOpt(ch,ff),Sig0_SOpt(ch,ff)]=MutualInfo_NMS(a(I_s & I_single ),s(I_s & I_single ),7,7);
       [Inf_RS_COpt(ch,ff),Sig0_COpt(ch,ff)]=MutualInfo_NMS(a(I_s & I_choice),s(I_s & I_choice),7,7);

   end
end  


figure(1)
for ch=1:size(Spectrum.TargetFFT)
subplot(4,4,ch)
plot(f_d(f_d<200),Inf_RS_CD(ch,:));hold on;
title('C Dir')
axis([0 200 0 0.2]);

subplot(4,4,ch+4)
plot(f_d(f_d<200),Inf_RS_COpt(ch,:));hold on;
title('C Opt')
axis([0 200 0 0.2]);
end

SpectrumInfoTarg.Inf_RS_SD=Inf_RS_SD;
SpectrumInfoTarg.Inf_RS_SOpt=Inf_RS_SOpt;
SpectrumInfoTarg.Inf_RS_CD=Inf_RS_CD;
SpectrumInfoTarg.Inf_RS_COpt=Inf_RS_COpt;
SpectrumInfoTarg.Sig_SD=Sig0_SD;
SpectrumInfoTarg.Sig_SOpt=Sig0_SOpt;
SpectrumInfoTarg.Sig_CD=Sig0_CD;
SpectrumInfoTarg.Sig_COpt=Sig0_COpt;
%%
%%%%%%%%%%%%%%%%%%%%%%%% Result %%%%%%%%%%%%%%%%%%%%%%%%%%
clear A A0 a
A=Spectrum.ResultFFT_exp; A=permute(A,[  2 1 3]);
A0=downsample(A,2);
A0=permute(A0,[2 1 3]);
A0=A0(:,f_d<200,:);



clear Inf_RS*
%%%%  For reward prediction
s=table2array(TrialInfo(:,14));I_s=(s>0);
%ss=(s>1)+2*(s==1);
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_Rexp(ch,ff),Sig0_Rexp(ch,ff)]=MutualInfo_NMS(a(I_s),s(I_s),4,4);
   end
end


clear A A0 a
A=Spectrum.ResultFFT; A=permute(A,[  2 1 3]);
A0=downsample(A,2);
A0=permute(A0,[2 1 3]);
A0=A0(:,f_d<200,:);


%%%%  For result   Reward amount
s=table2array(TrialInfo(:,17));I_s=(s>0);
%ss=(s>1)+2*(s==1);
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_R(ch,ff),Sig0_R(ch,ff)]=MutualInfo_NMS(a(I_s),s(I_s),4,4);
   end
end


%%%%  For result   Win/Loss
s=table2array(TrialInfo(:,17));I_s=(s>0);
ss=(s>1)+2*(s==1);
%%%%% For win/loss
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_W(ch,ff),Sig0_W(ch,ff)]=MutualInfo_NMS(a(I_s ),ss(I_s ),2,2);
   end
end


for ch=1:size(Spectrum.TargetFFT)
subplot(4,4,ch+8)
plot(f_d(f_d<200),Inf_RS_Rexp(ch,:));hold on;
title('RewardExp')
axis([0 200 0 0.2]);

subplot(4,4,ch+12)
plot(f_d(f_d<200),Inf_RS_R(ch,:));hold on;
title('Result')
axis([0 200 0 0.2]);
end


SpectrumInfoResult.Inf_RS_W=Inf_RS_W;
SpectrumInfoResult.Inf_RS_R=Inf_RS_R;
SpectrumInfoResult.Inf_RS_Rexp=Inf_RS_Rexp;

SpectrumInfoResult.Sig_W=Sig0_W;
SpectrumInfoResult.Sig_R=Sig0_R;
SpectrumInfoResult.Sig_Rexp=Sig0_Rexp;

% d=strfind(name,'.')-1;
h=figure(1);
print( h, '-djpeg', [Dirfigure,name,'MI']);
% h=figure(2)
% print( h, '-djpeg', [Dirfigure,name(1:d),'MPCh2']);
% SpectrumInfo.Inf_RS_norm_D=Inf_RS_norm_D;
% SpectrumInfo.Inf_RS_inact_D=Inf_RS_inact_D;
save(file,'SpectrumInfoResult','SpectrumInfoTarg','-append');