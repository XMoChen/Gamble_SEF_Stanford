
clear all
Dir='D:\Projects\SEFCooling\data\';
fileDir=dir(Dir);
fileDir_name={fileDir(4:25).name};
Dirfigure='D:\Projects\SEFCooling\figures\';

close all;
i=3;%[1:length(fileDir_name)]  %  9
file=[Dir,fileDir_name{i}];


load(file,'Infortable')
I_choice=table2array(Infortable(:,3))~=38 & table2array(Infortable(:,8))~=0 & table2array(Infortable(:,9))~=0 & table2array(Infortable(:,6))>0;
Infortable0=table2array(Infortable);
table_choice=Infortable0(I_choice,:);
I_control=table_choice(:,1)>35;
I_inact=table_choice(:,1)<15;


xdata(1,:)=[0.4 0.8 0.2 0.4 0.8 0.2 0.4 ];
xdata(2,:)=[3 3 5 5 5 9 9];
Prob=xdata(1,:);%[0.2 1 0.8 0.2 0.2 1 0.8 0.2 0.4];
Vmax=xdata(2,:);%[3 3 5 5 5 9 9];
V=Prob.*Vmax+(1-Prob)*1;
CV=round(sqrt((Prob.*(Vmax-V).^2+(1-Prob).*(1-V).^2))./V*100)/100;
Var=round(sqrt((Prob.*(Vmax-V).^2+(1-Prob).*(1-V).^2))*100)/100;



for tr=1:size(table_choice,1)
    
    Opt1_V=V(table_choice(tr,7));
    Opt2_V=V(table_choice(tr,8));
    Opt1_CV=Var(table_choice(tr,7));
    Opt2_CV=Var(table_choice(tr,8));
    
    if (Opt1_V>=Opt2_V & table_choice(tr,6)==1) | (Opt1_V<=Opt2_V & table_choice(tr,6)==2)
        HighV_I(tr)=1;
    elseif (Opt1_V>Opt2_V & table_choice(tr,6)==2) | (Opt1_V<Opt2_V & table_choice(tr,6)==1)
        HighV_I(tr)=-1;
    end
    
    if (Opt1_CV>=Opt2_CV & table_choice(tr,6)==1) | (Opt1_CV<=Opt2_CV & table_choice(tr,6)==2)
        HighCV_I(tr)=1;
    elseif (Opt1_CV>Opt2_CV & table_choice(tr,6)==2) | (Opt1_CV<Opt2_CV & table_choice(tr,6)==1)
        HighCV_I(tr)=-1;
    end
end

% trialN=50;
% kernel=ones(trialN,1)'/50;
sigma=30;
kernel=HalfGaussian(sigma);
HighV_Index=HighV_I==1;
HighV_Perc=convn(HighV_Index',kernel','same')';
HighV_Perc_con=HighV_Perc;
HighV_Perc_con(~I_control)=nan;

% LowV=HighV_I==-1;
% sigma=10;
% kernel=HalfGaussian(sigma);
% HighV_Index_Ins=convn(HighV_Index',kernel','same')';
% LowV_Ins=convn(LowV',kernel,'same')';
% HighV_Perc=HighV_Index_Ins./(LowV_Ins);
% % LowV_Perc=LowV_Ins./(HighV_Index_Ins);
%
% % HighVPerc=zeros(size(table_choice,1),1);
% HighVPerc=atan(HighV_Perc')*180/pi;
subplot(221)
plot(HighV_Perc(100:1100),'Color',[1, 0.5,0]); hold on;
plot(HighV_Perc_con(100:1100),'Color',[0, 0.5,0]); hold on;
axis([1 1000 0 1]);
box off; set(gca,'TickDir','out');

HighCV_Index=HighCV_I==1;
% trialN=30;
% kernel=ones(1,trialN);
% LowV=HighV_I==-1;
%  sigma=50;
%  kernel=HalfGaussian(sigma);
subplot(223)
HighCV_Perc=convn(HighCV_Index',kernel','same')';
HighCV_Perc_con=HighCV_Perc;
HighCV_Perc_con(~I_control)=nan;


plot(HighCV_Perc(100:1100),'Color',[1, 0.5,0]);hold on;
plot(HighCV_Perc_con(100:1100),'Color',[0, 0.5,0]); hold on;
axis([1 1000 0 1]);
box off; set(gca,'TickDir','out');

