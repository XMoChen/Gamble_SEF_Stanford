function SpectMI(Dirfigure,name,file)
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
I_norm=table2array(Infortable(:,1))>35 & table2array(Infortable(:,3))<38;
I_cool=table2array(Infortable(:,1))<10 & table2array(Infortable(:,3))<38;

%%%%% For direction
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_norm_SD(ch,ff),Sig0_norm_SD(ch,ff)]=MutualInfo_NMS(a(I_s & I_norm & I_single),s(I_s & I_norm & I_single),4,4);
       [Inf_RS_inact_SD(ch,ff),Sig0_inact_SD(ch,ff)]=MutualInfo_NMS(a(I_s & I_cool & I_single),s(I_s & I_cool & I_single),4,4);

       [Inf_RS_norm_CD(ch,ff),Sig0_norm_CD(ch,ff)]=MutualInfo_NMS(a(I_s & I_norm & I_choice),s(I_s & I_norm & I_choice),4,4);
       [Inf_RS_inact_CD(ch,ff),Sig0_inact_CD(ch,ff)]=MutualInfo_NMS(a(I_s & I_cool & I_choice),s(I_s & I_cool & I_choice),4,4);


   end
end


%%%%% For Option
s=table2array(Infortable(:,9));I_s=(s>0);
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_norm_SOpt(ch,ff),Sig0_norm_SOpt(ch,ff)]=MutualInfo_NMS(a(I_s & I_norm & I_single),s(I_s & I_norm & I_single),4,4);
       [Inf_RS_inact_SOpt(ch,ff),Sig0_inact_SOpt(ch,ff)]=MutualInfo_NMS(a(I_s & I_cool & I_single),s(I_s & I_cool & I_single),4,4);

       [Inf_RS_norm_COpt(ch,ff),Sig0_norm_COpt(ch,ff)]=MutualInfo_NMS(a(I_s & I_norm & I_choice),s(I_s & I_norm & I_choice),4,4);
       [Inf_RS_inact_COpt(ch,ff),Sig0_inact_COpt(ch,ff)]=MutualInfo_NMS(a(I_s & I_cool & I_choice),s(I_s & I_cool & I_choice),4,4);


   end
end

SpectrumInfoTarg.Inf_RS_norm_SD=Inf_RS_norm_SD;
SpectrumInfoTarg.Inf_RS_inact_SD=Inf_RS_inact_SD;
SpectrumInfoTarg.Inf_RS_norm_CD=Inf_RS_norm_CD;
SpectrumInfoTarg.Inf_RS_inact_CD=Inf_RS_inact_CD;
SpectrumInfoTarg.Sig_norm_SD=Inf_RS_norm_SD;
SpectrumInfoTarg.Sig_inact_SD=Inf_RS_inact_SD;
SpectrumInfoTarg.Sig_norm_CD=Inf_RS_norm_CD;
SpectrumInfoTarg.Sig_inact_CD=Inf_RS_inact_CD;

% SpectrumInfoTarg.Inf_RS_norm_Opt=Inf_RS_norm_Opt;
% SpectrumInfoTarg.Inf_RS_inact_Opt=Inf_RS_inact_Opt;
SpectrumInfoTarg.Inf_RS_norm_SOpt=Inf_RS_norm_SOpt;
SpectrumInfoTarg.Inf_RS_inact_SOpt=Inf_RS_inact_SOpt;
SpectrumInfoTarg.Inf_RS_norm_COpt=Inf_RS_norm_COpt;
SpectrumInfoTarg.Inf_RS_inact_COpt=Inf_RS_inact_COpt;
SpectrumInfoTarg.Sig_norm_SOpt=Sig0_norm_SOpt;
SpectrumInfoTarg.Sig_inact_SOpt=Sig0_inact_SOpt;
SpectrumInfoTarg.Sig_norm_COpt=Sig0_norm_COpt;
SpectrumInfoTarg.Sig_inact_COpt=Sig0_inact_COpt;

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
%%
figure(1)
subplot(241)
plot(f_d(f_d<200),Inf_RS_norm_COpt(1,:));hold on;
plot(f_d(f_d<200),Inf_RS_inact_SOpt(1,:));hold on;
title('ch1 SD')
axis([0 200 0 0.2]);

subplot(242)
plot(f_d(f_d<200),Inf_RS_norm_SD(2,:));hold on;
plot(f_d(f_d<200),Inf_RS_inact_SD(2,:));hold on;
axis([0 200 0 0.2]);
title('ch2 SD')

subplot(243)
plot(f_d(f_d<200),Inf_RS_norm_CD(1,:));hold on;
plot(f_d(f_d<200),Inf_RS_inact_CD(1,:));hold on;
title('ch1 CD')
axis([0 200 0 0.2]);

subplot(244)
plot(f_d(f_d<200),Inf_RS_norm_CD(2,:));hold on;
plot(f_d(f_d<200),Inf_RS_inact_CD(2,:));hold on;  
axis([0 200 0 0.2]);
title('ch2 CD')
%%
%%%%%%%%%%%%%%%%%%%%%%%% Result %%%%%%%%%%%%%%%%%%%%%%%%%%
%%reward expectation
%%%%%
clear A A0 a
A=Spectrum.ResultExpFFT; A=permute(A,[  2 1 3]);
A0=downsample(A,2);
A0=permute(A0,[2 1 3]);
A0=A0(:,f_d<200,:);


% clear Inf_RS*
%%%%  For result   Option/Option
s=table2array(Infortable(:,9));I_s=(s>0);
%ss=(s>1)+2*(s==1);
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_norm_Exp(ch,ff),Sig0_norm_Exp(ch,ff)]=MutualInfo_NMS(a(I_s & I_norm),s(I_s & I_norm),4,4);
       [Inf_RS_inact_Exp(ch,ff),Sig0_inact_Exp(ch,ff)]=MutualInfo_NMS(a(I_s & I_cool),s(I_s & I_cool),4,4);

   end
end


%%%%%
clear A A0 a
A=Spectrum.ResultFFT; A=permute(A,[  2 1 3]);
A0=downsample(A,2);
A0=permute(A0,[2 1 3]);
A0=A0(:,f_d<200,:);



%%%%  For result   Option/Option
s=table2array(Infortable(:,20));I_s=(s>0);
%ss=(s>1)+2*(s==1);
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_norm_R(ch,ff),Sig0_norm_R(ch,ff)]=MutualInfo_NMS(a(I_s & I_norm),s(I_s & I_norm),4,4);
       [Inf_RS_inact_R(ch,ff),Sig0_inact_R(ch,ff)]=MutualInfo_NMS(a(I_s & I_cool),s(I_s & I_cool),4,4);

   end
end


%%%%  For result   Win/Loss
s=table2array(Infortable(:,20));I_s=(s>0);
ss=(s>1)+2*(s==1);
%%%%% For win/loss
for ch=1:size(Spectrum.TargetFFT)
   A0_c=squeeze( A0(ch,:,:));
   for ff=1:size(A0_c,1)
       a=A0_c(ff,:);
       [Inf_RS_norm_W(ch,ff),Sig0_norm_W(ch,ff)]=MutualInfo_NMS(a(I_s & I_norm),ss(I_s & I_norm),2,2);
       [Inf_RS_inact_W(ch,ff),Sig0_inact_W(ch,ff)]=MutualInfo_NMS(a(I_s & I_cool),ss(I_s & I_cool),2,2);

   end
end
%%
subplot(245)
plot(f_d(f_d<200),Inf_RS_norm_Exp(1,:));hold on;
plot(f_d(f_d<200),Inf_RS_inact_Exp(1,:));hold on;
legend({'N','In'});
axis([0 200 0 0.1]);
title('ch1 RewExp')

subplot(246)
plot(f_d(f_d<200),Inf_RS_norm_Exp(2,:));hold on;
plot(f_d(f_d<200),Inf_RS_inact_Exp(2,:));hold on;
axis([0 200 0 0.1]);
title('ch2 RewExp')

subplot(247)
plot(f_d(f_d<200),Inf_RS_norm_R(1,:));hold on;
plot(f_d(f_d<200),Inf_RS_inact_R(1,:));hold on;
axis([0 200 0 0.1]);
title('ch1 Reward')

subplot(248)
plot(f_d(f_d<200),Inf_RS_norm_R(2,:));hold on;
plot(f_d(f_d<200),Inf_RS_inact_R(2,:));hold on;
axis([0 200 0 0.1]);
title('ch2 Reward')


SpectrumInfoResult.Inf_RS_norm_W=Inf_RS_norm_W;
SpectrumInfoResult.Inf_RS_inact_W=Inf_RS_inact_W;
SpectrumInfoResult.Inf_RS_norm_R=Inf_RS_norm_R;
SpectrumInfoResult.Inf_RS_inact_R=Inf_RS_inact_R;
SpectrumInfoResult.Inf_RS_norm_Exp=Inf_RS_norm_Exp;
SpectrumInfoResult.Inf_RS_inact_Exp=Inf_RS_inact_Exp;
SpectrumInfoResult.Sig_norm_W=Sig0_norm_W;
SpectrumInfoResult.Sig_inact_W=Sig0_inact_W;
SpectrumInfoResult.Sig_norm_R=Sig0_norm_R;
SpectrumInfoResult.Sig_inact_R=Sig0_inact_R;
SpectrumInfoResult.Sig_norm_Exp=Sig0_norm_Exp;
SpectrumInfoResult.Sig_inact_Exp=Sig0_inact_Exp;

% d=strfind(name,'.')-1;
h=figure(1)
print( h, '-djpeg', [Dirfigure,name,'MI']);
% h=figure(2)
% print( h, '-djpeg', [Dirfigure,name(1:d),'MPCh2']);
% SpectrumInfo.Inf_RS_norm_D=Inf_RS_norm_D;
% SpectrumInfo.Inf_RS_inact_D=Inf_RS_inact_D;
save(file,'SpectrumInfoResult','SpectrumInfoTarg','-append');