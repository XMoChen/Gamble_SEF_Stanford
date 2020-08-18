clear all

close all



% FlagName={'cgo03082012-01','cgo03112012-01','cgo03132012-01','cgo03152012-01','cgo03202012-02','cgo03252012-01','cgo03272012-01','cgo03292012-01','cgo04012012-01','cgo04032012-01',...
%            'cgonew04242012-01','cgonew05102012-01','cgonew05122012','cgonew05132012-01','cgonew05152012-01',...
% 'I10222012cool','Icool10262012-01','Igocooling11022012_mrg','Igocooling11022012c','Igocool11122012-01','Igocool11142012-01','Igocool11142012b-01','Igocool11182012-01','Igocool11182012b','Igocool12032012','Icgo12052012new','Icgo12052012new(b)','Icgo12102012new_mrg','Icgo12122012new'...
% };


FlagName={'cgo03082012-01','cgo03112012-01','cgo03132012-01','cgo03152012-01','cgo03202012-02','cgo03252012-01','cgo03272012-01','cgo03292012-01','cgo04012012-01','cgo04032012-01','cgonew05102012-01_sub1','cgonew05102012-01_sub2','cgonew04242012-01_sub1','cgonew04242012-01_sub2','cgonew05122012','Ago03082013Newb','I10222012cool','Icool10262012-01','Igocooling11022012_mrg','Igocooling11022012c','Igocool11122012-01','Igocool11142012-01','Igocool11142012b-01','Igocool11182012-01','Igocool11182012b','Igocool12032012','Icgo12052012new','Icgo12052012new(b)','Icgo12102012new_mrg','Icgo12122012new','Icgo12122012new(b)'};
%    FlagName={ 'I10222012cool','Icool10262012-01','Igocooling11022012_mrg','Igocooling11022012c','Igocool11122012-01','Igocool11142012-01','Igocool11142012b-01','Igocool11182012-01','Igocool11182012b','Igocool12032012','Icgo12052012new','Icgo12052012new(b)','Icgo12102012new_mrg','Icgo12122012new','Icgo12122012new(b)'};

%  FlagName={};


Break_sac_Perc=zeros(31,2);
for i=1:length(FlagName)%[1:9]%1:14  18
    
    %%% expected value for each options
    if i<=10 || (i>16 & i<=26)
        Prob=[0.4 0.8 0.2 0.4 0.8 0.2 0.4 ];
        Vmax=[3 3 5 5 5 9 9];
    else
        Prob=[0.6 0.8 0.2 0.6 0.8 0.2 0.6];
        Vmax=[3 3 5 5 5 9 9];
    end
    V=Prob.*Vmax+(1-Prob)*1;
    
    file=char(FlagName(i));
    
    
    
    clear Eye* BreakSac_all*  
    
    
    if exist(['D:\Projects\SEFCooling\data\data\',file,'\',file,'LFP.mat'])~=0
        load(['D:\Projects\SEFCooling\data\data\',file,'\',file,'LFP.mat'],'EyeX','EyeY');
    end
    
    load(['D:\Projects\SEFCooling\data\data\',file,'\',file,'.mat'],'table','eventtimes')
    
    Sac_RT_all=nan(size(table,1),1);
    %     Drop_t=
    
    if exist('EyeX')~=0
        
        
        
        % load('/Users/xiaomo/Documents/PhD/gamble exp/data/Icool10262012-01/Icool10262012-01LFP.mat','EyeX','EyeY');
        
        % load('/Users/xiaomo/Documents/PhD/gamble exp/data/Icool10262012-01/Icool10262012-01.mat','table')
        
        
        WindowSize=20;
        
        b_filter=(1/WindowSize)*ones(1,WindowSize);
        
        a_filter=1;
        
        
        
        TargOn=round(table(:,14));
        
        TargIn=round(table(:,15));
        
        %         ResultOn=round(table(:,20));
        
        fix1_x=[];fix2_x=[];fix1_y=[];fix2_y=[];
        
        V0_Sac_1=[];V0_Sac_2=[];
        T_Sac_2=[];T_Sac_1=[];
        
        X0_Sac_41=[]; Y0_Sac_41=[];
        X0_Rew_41=[]; Y0_Rew_41=[];
        X0_NRew_41=[]; Y0_NRew_41=[];
        X0_Sac_42=[]; Y0_Sac_42=[];
        X0_Rew_42=[]; Y0_Rew_42=[];
        X0_NRew_42=[]; Y0_NRew_42=[];
        BreakSac_all=zeros(size(table,1),1);

        
        
        for d=1:4
            
            for temp=1:2
                
                clear I EyeX_I EyeY_I TargOn_I TargIn_I  T_Sac X_Reward* X_Sac* X0_Sac  Y_Reward* Y_Sac* Y0_Sac X0_Reward Y0_Reward
                
                if temp==1
                    I=table(:,27)==d & table(:,18)~=0 & table(:,17)~=0 & table(:,35)>30  & sum(table(:,7:8),2)==0;
                else
                    I=table(:,27)==d & table(:,18)~=0 & table(:,17)~=0 & table(:,35)<15  & sum(table(:,7:8),2)==0;
                end
                
                
                
                
                
                EyeX_I=EyeX(I,:);
                
                EyeY_I=EyeY(I,:);
                
                
                
                
                
                TargOn_I=TargOn(I )';
                
                TargIn_I=TargIn(I )';
                
                ResultOn_I=eventtimes(I,24)';
                
                RewardI=table(I,21);
                %                 for t=1:sum(I==1)
                %                     if RewardI(t)~=0
                %                     t_reward=ResultOn_I(t);
                %                     elseif RewardI(t)==0
                %                     ResultOn_I(t)=t_reward;
                %                     end
                %                 end
                Break_sac=zeros(1,sum(I==1));
                for t=1:sum(I==1)
                    if   ResultOn_I(t)>0  &  ResultOn_I(t)<size(EyeX_I,2)
                        X_Reward(t,:)=filtfilt(b_filter,a_filter,EyeX_I(t,(ResultOn_I(t)-700):(ResultOn_I(t)+600))-nanmean(EyeX_I(t,(TargIn_I(t)-200):(TargIn_I(t)-100)),2));
                        
                        Y_Reward(t,:)=filtfilt(b_filter,a_filter,EyeY_I(t,(ResultOn_I(t)-700):(ResultOn_I(t)+600))-nanmean(EyeY_I(t,(TargIn_I(t)-200):(TargIn_I(t)-100)),2));
                        
                        XY_Reward(t,:)=sqrt(X_Reward(t,:).^2+X_Reward(t,:).^2);
                        
                        t00=700;
                        if RewardI(t)==0
                            t00=400;
                            while (XY_Reward(t,t00)<(nanmean(XY_Reward(t,500:600))+30) & XY_Reward(t,t00)>(nanmean(XY_Reward(t,500:600))-30)) & t00<900
                                
                                t00=t00+1;
                            end
                            if t00==900
                                t00=700;
                            end
                            
                            if temp==1
                                if max(diff(XY_Reward(t,(t00-300):(t00+100))))*55>20
                                    Break_sac(t)=1;  % break by saccade during control
                                else
                                    Break_sac(t)=-1;
                                end
                            else
                                if max(diff(XY_Reward(t,(t00-300):(t00+100))))*55>20
                                    Break_sac(t)=2; % break by saccade during inactivation
                                else
                                    Break_sac(t)=-2;
                                end
                            end
                        end
                        % t00
                        X0_Reward(t,:)=X_Reward(t,(t00-300):(t00+200));
                        Y0_Reward(t,:)=Y_Reward(t,(t00-300):(t00+200));
                        
                        
                    else
                        X0_Reward(t,:)=nan(1,501);
                        Y0_Reward(t,:)=nan(1,501);
                    end
                    
                    if TargIn_I(t)<10000
                        
                        X_Sac(t,:)=filtfilt(b_filter,a_filter,EyeX_I(t,(TargIn_I(t)-700):(TargIn_I(t)+2200))-nanmean(EyeX_I(t,(TargIn_I(t)-200):(TargIn_I(t)-150)),2));
                        
                        Y_Sac(t,:)=filtfilt(b_filter,a_filter,EyeY_I(t,(TargIn_I(t)-700):(TargIn_I(t)+2200))-nanmean(EyeY_I(t,(TargIn_I(t)-200):(TargIn_I(t)-150)),2));
                        
                        XY_Sac(t,:)=sqrt(X_Sac(t,:).^2+Y_Sac(t,:).^2);
                        
                        
                        t0=520;
                        
                        while XY_Sac(t,t0)<(nanmean(XY_Sac(t,500:520))+100) & t0<900
                            
                            t0=t0+1;
                        end
                        
                        
                        
                        
                        T_Sac(t)=t0+TargIn_I(t)-700-TargOn_I(t);
                        T_Result(t)= round(ResultOn_I(t)-TargIn_I(t)+700);
                        
                        
                        
                        
                        t1=size(XY_Sac,2);
                        
                        while XY_Sac(t,t1)>(nanmean(XY_Sac(t,(end-1900):end-2000))-30) & t1>500
                            t1=t1-1;
                            
                        end
                        
                        
                        T_Sac_end(t)=t1;
                        rt0=T_Result(t);
                        
                        
                        if t0~=800 & t1>500
                            
                            Sac_end_V(t)=10*1000/(t1-t0);
                            
                            
                            
                            X0_Sac(t,:)=X_Sac(t,(t0-500):(t0+900));
                            
                            Y0_Sac(t,:)=Y_Sac(t,(t0-500):(t0+900));
                            
                            XY0_Sac(t,:)=XY_Sac(t,(t0-500):(t0+900));
                            
                            V0_Sac(t)=max(diff(XY_Sac(t,(t0-100):(t0+100))));
                            
                            
                        else
                            X0_Sac(t,:)=nan;
                            Y0_Sac(t,:)=nan;
                            XY0_Sac(t,:)=nan;
                            V0_Sac(t)=nan;
                        end
                        
                        
                        
                        
                    end
                    
                end
                
                %                 %%%% total number of breaking saccade
                %                 Break_sac_Perc(d,1)=Break_sac_Perc(d,1)+sum(Break_sac==1);
                %                 Break_sac_Perc(d,1)=Break_sac_Perc(d,1)+sum(Break_sac==2);
                
                uu0=0;
                for u1=1:7
                    a0(u1)=nanmean(T_Sac(table(I,23)==0 & table(I,24)==u1));
                    
                    for u2=1:7
                        uu0=uu0+1;
                        if sum(table(I,23)~=0 & table(I,24)==u1 )>=1
                            b0(uu0)=nanmean(T_Sac(table(I,23)~=0 & table(I,24)==u1 & table(I,25)==u2));
                        else
                            b0(uu0)=nan;
                        end
                    end
                end
                
                
                Sac_RT(d,temp,1)=nanmean(a0);
                Sac_RT(d,temp,2)=nanmean(b0);
                
                Sac_RT_all(I,1)=T_Sac;
                BreakSac_all(I,1)=Break_sac;
                
                clear I0
                
                Itr= X0_Sac(:,end)>(nanmean(X0_Sac(:,600))-0.5*nanstd(X0_Sac(:,600))) &  X0_Sac(:,600)<(nanmean(X0_Sac(:,600))+0.5*nanstd(X0_Sac(:,600))) ...
                    & Y0_Sac(:,600)>(nanmean(Y0_Sac(:,600))-0.5*nanstd(Y0_Sac(:,600))) &  Y0_Sac(:,600)<(nanmean(Y0_Sac(:,600))+0.5*nanstd(Y0_Sac(:,end)));
                %     & X0_Sac(:,t0-100)>(nanmean(X0_Sac(:,t0-100))-0.5*nanstd(X0_Sac(:,t0-100))) &  X0_Sac(:,t0-100)<(nanmean(X0_Sac(:,t0-100))+0.5*nanstd(X0_Sac(:,t0-100)))...
                %     & Y0_Sac(:,t0-100)>(nanmean(Y0_Sac(:,t0-100))-0.5*nanstd(Y0_Sac(:,t0-100))) &  Y0_Sac(:,t0-100)<(nanmean(Y0_Sac(:,t0-100))+0.5*nanstd(Y0_Sac(:,t0-100)))...
                %     & X0_Sac(:,t0-400)>(nanmean(X0_Sac(:,t0-400))-0.5*nanstd(X0_Sac(:,t0-400))) &  X0_Sac(:,t0-400)<(nanmean(X0_Sac(:,t0-400))+0.5*nanstd(X0_Sac(:,t0-400)))...
                %     & Y0_Sac(:,t0-400)>(nanmean(Y0_Sac(:,t0-400))-0.5*nanstd(Y0_Sac(:,t0-400))) &  Y0_Sac(:,t0-400)<(nanmean(Y0_Sac(:,t0-400))+0.5*nanstd(Y0_Sac(:,t0-400)));
                % %
                %
                
                %if temp==1
                
                I1=table(I,18)>0 & table(I,21)>0 ;
                I0=table(I,21)==0 & table(I,18)>0;
                
                
                dx=(10*sqrt(2)/2)/abs(nanmean(X0_Sac(I0 & Itr,600)));
                
                dy=(10*sqrt(2)/2)/abs(nanmean(Y0_Sac(I0 & Itr,600)));
                
                dxy=(10)/abs(nanmean(XY0_Sac(I0 & Itr,600)));
                
                %end
                
                X0_Sac=X0_Sac*dx;
                
                Y0_Sac=Y0_Sac*dy;
                
                V0_Sac=V0_Sac*dxy*1000;
                
                X_reward=X0_Reward*dx;
                Y_reward=Y0_Reward*dy;
                
                if temp==1
                    V0_Sac_1(d)=nanmean(V0_Sac(I1 & Itr));
                    X0_Sac_41=[X0_Sac_41;X0_Sac(I1 & Itr,:)];
                    Y0_Sac_41=[Y0_Sac_41;Y0_Sac(I1 & Itr,:)];
                    X0_Rew_41=[X0_Rew_41;X_reward(I1,:)];
                    Y0_Rew_41=[Y0_Rew_41;Y_reward(I1,:)];
                    X0_NRew_41=[X0_NRew_41;X_reward(I0,:)];
                    Y0_NRew_41=[Y0_NRew_41;Y_reward(I0,:)];
                else
                    V0_Sac_2(d)=nanmean(V0_Sac(I1 & Itr));
                    X0_Sac_42=[X0_Sac_42;X0_Sac(I1 & Itr,:)];
                    Y0_Sac_42=[Y0_Sac_42;Y0_Sac(I1 & Itr,:)];
                    X0_Rew_42=[X0_Rew_42;X_reward(I1,:)];
                    Y0_Rew_42=[Y0_Rew_42;Y_reward(I1,:)];
                    X0_NRew_42=[X0_NRew_42;X_reward(I0,:)];
                    Y0_NRew_42=[Y0_NRew_42;Y_reward(I0,:)];
                end
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%% example for trajectories
                %
                %                                   figure(1)
                % %
                %                                   PlotTrajectory;  %recording 18
                
            end
            
            
        end
        %%
        colormap hot
        xbins=-30:1:30;
        ybins=-30:1:30;
        FixControlH(i,:,:)=NormDis(D2Histogram([X0_Sac_41(:,400),Y0_Sac_41(:,400)],xbins,ybins));
        SacControlH(i,:,:)=NormDis(D2Histogram([X0_Sac_41(:,600),Y0_Sac_41(:,600)],xbins,ybins)/length(X0_Sac_41(:,400)~=0));
        
        FixControlH1d(i,:)=NormDis(hist(sqrt(X0_Sac_41(:,400).^2+Y0_Sac_41(:,400).^2),xbins)/length(X0_Sac_41(:,400)~=0));
        SacControlH1d(i,:)=NormDis(hist(sqrt(X0_Sac_41(:,600).^2+Y0_Sac_41(:,600).^2),xbins)/length(X0_Sac_41(:,400)~=0));
        
        SucControlH(i,:,:)=D2Histogram([X0_Rew_41(:,200),Y0_Rew_41(:,200)],xbins,ybins)/(length(X0_Rew_41(:,200)~=0)+length(X0_NRew_41(:,200)~=0));
        SucControlH1d(i,:,:)=NormDis(hist(sqrt(X0_Rew_41(:,200).^2+Y0_Rew_41(:,200).^2),xbins)/(length(X0_Rew_41(:,200)~=0)));
        FaiControlH(i,:,:)=D2Histogram([X0_NRew_41(:,350),Y0_NRew_41(:,350)],xbins,ybins)/(length(X0_Rew_41(:,200)~=0)+length(X0_NRew_41(:,200)~=0));
        
        FixInactH(i,:,:)=NormDis(D2Histogram([X0_Sac_42(:,400),Y0_Sac_42(:,400)],xbins,ybins)/length(X0_Sac_42(:,400)~=0));
        SacInactH(i,:,:)=NormDis(D2Histogram([X0_Sac_42(:,600),Y0_Sac_42(:,600)],xbins,ybins)/length(X0_Sac_42(:,400)~=0));
        FixInactH1d(i,:)=NormDis(hist(sqrt(X0_Sac_42(:,400).^2+Y0_Sac_42(:,400).^2),ybins)/length(X0_Sac_42(:,400)~=0));
        SacInactH1d(i,:)=NormDis(hist(sqrt(X0_Sac_42(:,600).^2+Y0_Sac_42(:,600).^2),ybins)/length(X0_Sac_42(:,400)~=0));
        SucInactH(i,:,:)=D2Histogram([X0_Rew_42(:,200),Y0_Rew_42(:,200)],xbins,ybins)/(length(X0_Rew_42(:,200)~=0)+length(X0_NRew_42(:,200)~=0));
        
        SucInactH1d(i,:,:)=NormDis(hist(sqrt(X0_Rew_42(:,200).^2+Y0_Rew_42(:,200).^2),xbins)/(length(X0_Rew_42(:,200)~=0)));
        FaiInactH(i,:,:)=D2Histogram([X0_NRew_42(:,350),Y0_NRew_42(:,350)],xbins,ybins)/(length(X0_Rew_42(:,200)~=0)+length(X0_NRew_42(:,200)~=0));
        
        FixateControlDelta(i)=nanstd(sqrt(X0_Sac_41(:,400).^2+Y0_Sac_41(:,400).^2));%/sqrt(length(X0_Sac_41(:,400)~=0))-1);
        SacControlDelta(i)=nanstd(sqrt(X0_Sac_41(:,600).^2+Y0_Sac_41(:,600).^2));%/sqrt(length(X0_Sac_41(:,600)~=0))-1);
        SucControlDelta(i,:,:)=nanstd(sqrt(X0_Rew_41(:,200).^2+Y0_Rew_41(:,200).^2));%/sqrt((length(X0_Rew_41(:,200)~=0))-1);
        
        FixateInactDelta(i)=nanstd(sqrt(X0_Sac_42(:,400).^2+Y0_Sac_42(:,400).^2));%/sqrt(length(X0_Sac_42(:,400)~=0))-1);
        SacInactDelta(i)=nanstd(sqrt(X0_Sac_42(:,600).^2+Y0_Sac_42(:,600).^2));%/sqrt(length(X0_Sac_42(:,600)~=0))-1);
        SucInactDelta(i,:,:)=nanstd(sqrt(X0_Rew_42(:,200).^2+Y0_Rew_42(:,200).^2));%/sqrt((length(X0_Rew_42(:,200)~=0))-1);
        
        
        
        
        %%
        
        
        V_S(i,:)=[nanmean(V0_Sac_1),nanmean(V0_Sac_2)];
        T_all(i,:)=[nanmean(T_Sac_1),nanmean(T_Sac_2)];
        
        I0= table(:,18)~=0 & table(:,17)~=0 & (table(:,35)>35) & sum(table(:,7:8),2)==0 & table(:,18)>0;
        ResultBreak(i,1)=sum(BreakSac_all(I0)==1)/sum(I0);
        I1= table(:,18)~=0 & table(:,17)~=0 & (table(:,35)<15) & sum(table(:,7:8),2)==0 & table(:,18)>0;
        ResultBreak(i,2)=sum(BreakSac_all(I1)==2)/sum(I1);
        
        
        dif_fix_x=nanstd(fix1_x)-nanstd(fix2_x);
        
        dif_fix_y=nanstd(fix1_y)-nanstd(fix2_y);
        
        z_x=[fix1_x,fix2_x];
        
        z_y=[fix1_y,fix2_y];
        
        
        
        Dif_fix_x(i,:)=[nanstd(fix1_x),nanstd(fix2_x)];
        
        Dif_fix_y(i,:)=[nanstd(fix1_y),nanstd(fix2_y)];
        
        
        
    else
        V_S(i,:)=nan;
        ResultBreak(i,1)=nan;
        ResultBreak(i,2)=nan;
        
        FixControlH(i,:,:)=nan;
        SacControlH(i,:,:)=nan;
        FixControlH1d(i,:)=nan;
        SacControlH1d(i,:)=nan;
        SucControlH(i,:,:)=nan;
        FaiControlH(i,:,:)=nan;
        
        FixInactH(i,:,:)=nan;
        SacInactH(i,:,:)=nan;
        FixInactH1d(i,:)=nan;
        SacInactH1d(i,:)=nan;
        SucInactH(i,:,:)=nan;
        FaiInactH(i,:,:)=nan;
        
        
        
        T_Sac=table(:,15)-table(:,14);
        
        for d=1:4
            for temp=1:2
                if temp==1
                    I=table(:,27)==d & table(:,18)~=0 & table(:,17)~=0 & table(:,35)>35;
                else
                    I=table(:,27)==d & table(:,18)~=0 & table(:,17)~=0 & table(:,35)<15;
                end
                
                uu0=0;
                for u1=1:7
                    a0(u1)=nanmean(T_Sac(table(I,23)==0 & table(I,24)==u1));
                    
                    for u2=1:7
                        uu0=uu0+1;
                        if sum(table(I,23)~=0 & table(I,24)==u1 & table(I,25)==u2)>=1
                            b0(uu0)=nanmean(T_Sac(table(I,23)~=0 & table(I,24)==u1 & table(I,25)==u2));
                        else
                            b0(uu0)=nan;
                        end
                    end
                end
                Sac_RT(d,temp,1)=nanmean(a0);
                Sac_RT(d,temp,2)=nanmean(b0);
            end
        end
    end
    
    
    
    Sac_setion_noChoice(i,:,:)=Sac_RT(:,:,1);
    Sac_setion_Choice(i,:,:)=Sac_RT(:,:,2);
    
    
    %%%% calculate drop rate by saccade
    clear ExpDif
    for trial=1:size(table,1)
        if table(trial,19)>0 & table(trial,18)==1 & table(trial,24)~=0 & table(trial,24)<8
            ExpDif(trial)=Vmax(table(trial,24))-V(table(trial,24));
        elseif  table(trial,19)>0 & table(trial,18)==2  & table(trial,24)~=0 & table(trial,24)<8
            ExpDif(trial)=1-V(table(trial,24));
        else
            ExpDif(trial)=nan;
        end
    end
    
    ExpDif_I=[-4.80000000000000,-3.20000000000000,-2.40000000000000,-1.60000000000000,-1.60000000000000,-1.20000000000000,-0.800000000000000,-0.800000000000000,0,0.399999999999999,0.800000000000000,0.800000000000000,1.20000000000000,1.60000000000000,2.40000000000000,3.20000000000000,4.80000000000000,6.40000000000000];%unique(ExpDif(ExpDif>=-100));
    u=0;
    for Vexp_dif=1:length(ExpDif_I)
        u=u+1;
        ved=ExpDif_I(u);
        if exist('BreakSac_all')==0
            Drop_Result(i,u,1)=sum( table(:,19)>0 & table(:,20)==0 & table(:,35)>30 & ExpDif'==ved)/sum(table(:,19)>0 & table(:,35)>30 & ExpDif'==ved);
            Drop_Result(i,u,2)=sum(table(:,19)>0 & table(:,20)==0 & table(:,35)<15 & ExpDif'==ved)/sum(table(:,19)>0 & table(:,35)<15 & ExpDif'==ved);
            Drop_Result(i,u,3)=ved;
            
            % Choice
            Drop_ResultCNC(i,u,1)=sum( table(:,16)>0 & table(:,19)>0 & table(:,20)==0 & table(:,23)~=0 & table(:,35)>30 & ExpDif'==ved)/sum(table(:,16)>0 & table(:,19)>0 & table(:,23)~=0 & table(:,35)>30 & ExpDif'==ved);
            Drop_ResultCNC(i,u,2)=sum(table(:,16)>0 & table(:,19)>0 & table(:,20)==0  & table(:,23)~=0 & table(:,35)<15 & ExpDif'==ved)/sum(table(:,16)>0 & table(:,19)>0 & table(:,23)~=0 & table(:,35)<15 & ExpDif'==ved);
            Drop_ResultCNC(i,u,3)=ved;
            %NoChoice
            Drop_ResultCNC(i,u,4)=sum(table(:,16)>0 &  table(:,19)>0 & table(:,20)==0  & table(:,23)==0 & table(:,35)>30 & ExpDif'==ved)/sum(table(:,16)>0 & table(:,19)>0 & table(:,23)==0 & table(:,35)>30 & ExpDif'==ved);
            Drop_ResultCNC(i,u,5)=sum(table(:,16)>0 & table(:,19)>0 & table(:,20)==0 &  table(:,23)==0 & table(:,35)<15 & ExpDif'==ved)/sum(table(:,16)>0 & table(:,19)>0 & table(:,23)==0 & table(:,35)<15 & ExpDif'==ved);
            Drop_ResultCNC(i,u,6)=ved;

  drp(i,2)= sum(table(:,16)>0 &  table(:,19)>0 & table(:,20)==0  & table(:,23)~=0 & table(:,35)<15)/sum(table(:,16)>0 &  table(:,19)>0  & table(:,23)~=0 & table(:,35)<15);
  drp(i,1)=sum(table(:,16)>0 &  table(:,19)>0 & table(:,20)==0  & table(:,23)~=0 & table(:,35)>30)/sum(table(:,16)>0 &  table(:,19)>0  & table(:,23)~=0 & table(:,35)>30);
  drp(i,3)=sum(table(:,16)>0 &  table(:,19)>0 & table(:,20)==0  & table(:,23)==0 & table(:,35)>30)/sum(table(:,16)>0 &  table(:,19)>0  & table(:,23)==0 & table(:,35)>30);
  drp(i,4)=sum(table(:,16)>0 &  table(:,19)>0 & table(:,20)==0  & table(:,23)==0 & table(:,35)<15)/sum(table(:,16)>0 &  table(:,19)>0  & table(:,23)==0 & table(:,35)<15);

            
            
        else
            Drop_Result(i,u,1)=sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0 & table(:,20)==0 & table(:,35)>30 & ExpDif'==ved)/sum(table(:,19)>0 & table(:,35)>30 & ExpDif'==ved);
            Drop_Result(i,u,2)=sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0 & table(:,20)==0 & table(:,35)<15 & ExpDif'==ved)/sum(table(:,19)>0 & table(:,35)<15 & ExpDif'==ved);
            Drop_Result(i,u,3)=ved;
            
            % Choice
            Drop_ResultCNC(i,u,1)=sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0  & table(:,23)~=0 & table(:,20)==0 & table(:,35)>30 & ExpDif'==ved)/sum(table(:,16)>0 & table(:,19)>0  & table(:,23)~=0 & table(:,35)>30 & ExpDif'==ved);
            Drop_ResultCNC(i,u,2)=sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0  & table(:,23)~=0 & table(:,20)==0 & table(:,35)<15 & ExpDif'==ved)/sum(table(:,16)>0 & table(:,19)>0 & table(:,23)~=0 & table(:,35)<15 & ExpDif'==ved);
            Drop_ResultCNC(i,u,3)=ved;
             
            % NoChoice           
            Drop_ResultCNC(i,u,4)=sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0 & table(:,20)==0  & table(:,23)==0 & table(:,35)>30 & ExpDif'==ved)/sum(table(:,16)>0 & table(:,19)>0 & table(:,23)==0 & table(:,35)>30 & ExpDif'==ved);
            Drop_ResultCNC(i,u,5)=sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0 & table(:,20)==0  & table(:,23)==0 & table(:,35)<15 & ExpDif'==ved)/sum(table(:,16)>0 & table(:,19)>0 & table(:,23)==0 & table(:,35)<15 & ExpDif'==ved);
            Drop_ResultCNC(i,u,6)=ved;

  drp(i,2)= sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0 & table(:,20)==0  & table(:,23)~=0 & table(:,35)<15)/sum( table(:,16)>0 & table(:,19)>0  & table(:,23)~=0 & table(:,35)<15);
  drp(i,1)=sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0 & table(:,20)==0  & table(:,23)~=0 & table(:,35)>30)/sum( table(:,16)>0 & table(:,19)>0  & table(:,23)~=0 & table(:,35)>30);
  drp(i,3)=sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0 & table(:,20)==0  & table(:,23)==0 & table(:,35)>30)/sum( table(:,16)>0 & table(:,19)>0  & table(:,23)==0 & table(:,35)>30);
  drp(i,4)=sum(table(:,16)>0 & BreakSac_all>0 & table(:,19)>0 & table(:,20)==0  & table(:,23)==0 & table(:,35)<15)/sum(table(:,16)>0 &  table(:,19)>0  & table(:,23)==0 & table(:,35)<15);

        end
    end
    
  
     save(['D:\Projects\SEFCooling\data\data\',file,'\',file,'.mat'],'Sac_RT_all','Sac_RT','-append');
end
save('sum2.mat');

%% Eeror trials by saccades
load('sum2.mat')
figure(7)
subplot(221)
%plot(squeeze(Drop_Result_g(1,:,3)),squeeze(nanmean(Drop_Result_g(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_Result(1,:,3)),squeeze(nanmean(Drop_Result(1:10,:,1),1)),squeeze(nanstd(Drop_Result(1:10,:,1),1)/sqrt(10-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
%plot(squeeze(Drop_Result_n(1,:,3)),squeeze(nanmean(Drop_Result_n(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_Result(11,:,3)),squeeze(nanmean(Drop_Result(11:16,:,1),1)),squeeze(nanstd(Drop_Result(11:16,:,1),1)/sqrt(6-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_n=[squeeze(nanmean(Drop_Result(1:10,:,1),1)),squeeze(nanmean(Drop_Result(11:16,:,1),1))];
drop_all_expdif=[squeeze(nanmean(Drop_Result(1:10,:,3),1)),squeeze(nanmean(Drop_Result(11:16,:,3),1))];

errorbar(squeeze(Drop_Result(1,:,3)),squeeze(nanmean(Drop_Result(1:10,:,2),1)),squeeze(nanstd(Drop_Result(1:10,:,2),1)/sqrt(10-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
errorbar(squeeze(Drop_Result(11,:,3)),squeeze(nanmean(Drop_Result(11:16,:,2),1)),squeeze(nanstd(Drop_Result(11:16,:,2),1)/sqrt(6-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_i=[squeeze(nanmean(Drop_Result(1:10,:,2),1)),squeeze(nanmean(Drop_Result(11:16,:,2),1))];
drop_all_expdif=[squeeze(nanmean(Drop_Result(1:10,:,3),1)),squeeze(nanmean(Drop_Result(11:16,:,3),1))];

exp_def = fittype('(a*exp(-a*(x)))');
II=~isnan(drop_all_n);
f_1 = fit(drop_all_expdif(II)',drop_all_n(II)',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1, 0]);%,'Lower', [0, -inf, 0, -inf]);
[~,aa_n]=predint(f_1,-5:0.1:7);
coeffvalues(f_1)
plot(-5:0.1:7,aa_n,'Color',[0, 0.5,0],'linewidth',1);hold on;
f_2 = fit(drop_all_expdif(~isnan(drop_all_i))',drop_all_i(~isnan(drop_all_i))',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1,  0]);%,,'Lower', [0, -inf, 0, -inf]);
[~,aa_i]=predint(f_2,-5:0.1:7);
coeffvalues(f_2)
plot(-5:0.1:7,aa_i,'Color',[1, 0.5,0],'linewidth',1)
axis([-6 8 -0.1 1])
box off;set(gca,'TickDir','out')
axis square


subplot(222)
%plot(squeeze(Drop_Result_g(1,:,3)),squeeze(nanmean(Drop_Result_g(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_Result(1,:,3)),squeeze(nanmean(Drop_Result(17:26,:,1),1)),squeeze(nanstd(Drop_Result(17:26,:,1),1)/sqrt(10-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
%plot(squeeze(Drop_Result_n(1,:,3)),squeeze(nanmean(Drop_Result_n(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_Result(11,:,3)),squeeze(nanmean(Drop_Result(27:end,:,1),1)),squeeze(nanstd(Drop_Result(27:end,:,1),1)/sqrt(6-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_n=[squeeze(nanmean(Drop_Result(17:26,:,1),1)),squeeze(nanmean(Drop_Result(27:end,:,1),1))];
drop_all_expdif=[squeeze(nanmean(Drop_Result(17:26,:,3),1)),squeeze(nanmean(Drop_Result(27:end,:,3),1))];

errorbar(squeeze(Drop_Result(1,:,3)),squeeze(nanmean(Drop_Result(17:26,:,2),1)),squeeze(nanstd(Drop_Result(17:26,:,2),1)/sqrt(10-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
errorbar(squeeze(Drop_Result(11,:,3)),squeeze(nanmean(Drop_Result(27:end,:,2),1)),squeeze(nanstd(Drop_Result(27:end,:,2),1)/sqrt(6-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_i=[squeeze(nanmean(Drop_Result(17:26,:,2),1)),squeeze(nanmean(Drop_Result(27:end,:,2),1))];
drop_all_expdif=[squeeze(nanmean(Drop_Result(17:26,:,3),1)),squeeze(nanmean(Drop_Result(27:end,:,3),1))];

exp_def = fittype('(a*exp(-a*(x)))');
II=~isnan(drop_all_n);
f_1 = fit(drop_all_expdif(II)',drop_all_n(II)',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1, 0]);%,'Lower', [0, -inf, 0, -inf]);
[~,aa_n]=predint(f_1,-5:0.1:7);
coeffvalues(f_1)
plot(-5:0.1:7,aa_n,'Color',[0, 0.5,0],'linewidth',1);hold on;
f_2 = fit(drop_all_expdif(~isnan(drop_all_i))',drop_all_i(~isnan(drop_all_i))',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1,  0]);%,,'Lower', [0, -inf, 0, -inf]);
[~,aa_i]=predint(f_2,-5:0.1:7);
coeffvalues(f_2)
plot(-5:0.1:7,aa_i,'Color',[1, 0.5,0],'linewidth',1)
axis([-6 8 -0.1 0.8])
box off;set(gca,'TickDir','out')
axis square




exp_def = fittype('(a*exp(-a*(x)))');

for ii=1:size(Drop_Result)
    a=Drop_Result(ii,:,3);
    b=Drop_Result(ii,:,2);
    c=Drop_Result(ii,:,1);
    II=a~=0 & ~isnan(b) & ~isnan(c);
    f1 = fit(squeeze(Drop_Result(ii,II,3))',squeeze(Drop_Result(ii,II,1))',exp_def,'Robust','LAR');%,'Upper', [10, -8],'Lower', [0.5, -20]);
    f2 = fit(squeeze(Drop_Result(ii,II,3))',squeeze(Drop_Result(ii,II,2))',exp_def,'Robust','LAR');%,'Upper', [10, -8],'Lower', [0.5, -20]);

    Drop_Result_con_coef(ii,:)=coeffvalues(f1);
    Drop_Result_inact_coef(ii,:)=coeffvalues(f2);

end


% Drop_Result_con_coef(Drop_Result_con_coef<-24.5)=nan;
% Drop_Result_inact_coef(Drop_Result_inact_coef<-24.5)=nan;

subplot(223)
a=[nanmean(Drop_Result_con_coef(1:16,1)),nanmean(Drop_Result_inact_coef(1:16,1));...
    nanmean(Drop_Result_con_coef(17:end,1)),nanmean(Drop_Result_inact_coef(17:end,1))];
a_nanstd=[nanstd(Drop_Result_con_coef(1:16,1))/sqrt(16-1),nanstd(Drop_Result_inact_coef(1:16,1))/sqrt(16-1);...
    nanstd(Drop_Result_con_coef(17:end,1))/sqrt(15-1),nanstd(Drop_Result_inact_coef(17:end,1))/sqrt(15-1)];
errorbarplot(a,a_nanstd);
box off;set(gca,'TickDir','out')
axis([0.5 2.5 0 0.25]);
[~,p(3)]=ttest(Drop_Result_con_coef,Drop_Result_inact_coef)
[~,p(1)]=ttest(Drop_Result_con_coef(1:16,1),Drop_Result_inact_coef(1:16,1));
[~,p(2)]=ttest(Drop_Result_con_coef(17:end,1),Drop_Result_inact_coef(17:end,1));
axis square


subplot(224)
plot((ones(16,1)*[1,2])',([Drop_Result_con_coef(1:16,1),Drop_Result_inact_coef(1:16,1)])','k','MarkerSize',15);hold on
plot((ones(15,1)*[3,4])',([Drop_Result_con_coef(17:end,1),Drop_Result_inact_coef(17:end,1)])','k','MarkerSize',15);hold on
plot(ones(1,16),Drop_Result_con_coef(1:16,1),'g.','MarkerSize',15);hold on
plot(2*ones(1,16),Drop_Result_inact_coef(1:16,1),'b.','MarkerSize',15);hold on
plot(3*ones(1,15),Drop_Result_con_coef(17:end,1),'g.','MarkerSize',15);hold on
plot(4*ones(1,15),Drop_Result_inact_coef(17:end,1),'b.','MarkerSize',15);hold on


[~,p(3)]=ttest(Drop_Result_con_coef,Drop_Result_inact_coef)
[~,p(1)]=ttest(Drop_Result_con_coef(1:16,1),Drop_Result_inact_coef(1:16,1));
[~,p(2)]=ttest(Drop_Result_con_coef(17:end,1),Drop_Result_inact_coef(17:end,1));
axis square
box off;set(gca,'TickDir','out')
axis([0 5 0 0.3]);
%figure()
% plot(ones(1,16),Drop_Result_con_coef(1:16,1),'g.','MarkerSize',15);hold on
% plot(2*ones(1,16),Drop_Result_inact_coef(1:16,1),'b.','MarkerSize',15);hold on
% plot(3*ones(1,15),Drop_Result_con_coef(17:end,1),'g.','MarkerSize',15);hold on
% plot(4*ones(1,15),Drop_Result_inact_coef(17:end,1),'b.','MarkerSize',15);hold on
% 
% [~,p(3)]=ttest(Drop_Result_con_coef,Drop_Result_inact_coef)
% [~,p(1)]=ttest(Drop_Result_con_coef(1:16,1),Drop_Result_inact_coef(1:16,1));
% [~,p(2)]=ttest(Drop_Result_con_coef(17:end,1),Drop_Result_inact_coef(17:end,1));
% axis square

%%

%%%%% error rate seperate by 
figure(8)
subplot(241)
%plot(squeeze(Drop_Result_g(1,:,3)),squeeze(nanmean(Drop_Result_g(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_ResultCNC(1,:,3)),squeeze(nanmean(Drop_ResultCNC(1:10,:,1),1)),squeeze(nanstd(Drop_ResultCNC(1:10,:,1),1)/sqrt(10-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
%plot(squeeze(Drop_ResultCNC_n(1,:,3)),squeeze(nanmean(Drop_ResultCNC_n(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_ResultCNC(11,:,3)),squeeze(nanmean(Drop_ResultCNC(11:16,:,1),1)),squeeze(nanstd(Drop_ResultCNC(11:16,:,1),1)/sqrt(6-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_n=[squeeze(nanmean(Drop_ResultCNC(1:10,:,1),1)),squeeze(nanmean(Drop_ResultCNC(11:16,:,1),1))];
drop_all_expdif=[squeeze(nanmean(Drop_ResultCNC(1:10,:,3),1)),squeeze(nanmean(Drop_ResultCNC(11:16,:,3),1))];

errorbar(squeeze(Drop_ResultCNC(1,:,3)),squeeze(nanmean(Drop_ResultCNC(1:10,:,2),1)),squeeze(nanstd(Drop_ResultCNC(1:10,:,2),1)/sqrt(10-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
errorbar(squeeze(Drop_ResultCNC(11,:,3)),squeeze(nanmean(Drop_ResultCNC(11:16,:,2),1)),squeeze(nanstd(Drop_ResultCNC(11:16,:,2),1)/sqrt(6-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_i=[squeeze(nanmean(Drop_ResultCNC(1:10,:,2),1)),squeeze(nanmean(Drop_ResultCNC(11:16,:,2),1))];
drop_all_expdif=[squeeze(nanmean(Drop_ResultCNC(1:10,:,3),1)),squeeze(nanmean(Drop_ResultCNC(11:16,:,3),1))];

exp_def = fittype('(a*exp(-a*(x)))');
II=~isnan(drop_all_n);
f_1 = fit(drop_all_expdif(II)',drop_all_n(II)',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1, 0]);%,'Lower', [0, -inf, 0, -inf]);
[~,aa_n]=predint(f_1,-5:0.1:7);
coeffvalues(f_1)
plot(-5:0.1:7,aa_n,'Color',[0, 0.5,0],'linewidth',1);hold on;
f_2 = fit(drop_all_expdif(~isnan(drop_all_i))',drop_all_i(~isnan(drop_all_i))',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1,  0]);%,,'Lower', [0, -inf, 0, -inf]);
[~,aa_i]=predint(f_2,-5:0.1:7);
coeffvalues(f_2)
plot(-5:0.1:7,aa_i,'Color',[1, 0.5,0],'linewidth',1)
axis([-6 8 -0.1 1])
box off;set(gca,'TickDir','out')
axis square


subplot(242)
%plot(squeeze(Drop_ResultCNC_g(1,:,3)),squeeze(nanmean(Drop_ResultCNC_g(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_ResultCNC(1,:,3)),squeeze(nanmean(Drop_ResultCNC(17:26,:,1),1)),squeeze(nanstd(Drop_ResultCNC(17:26,:,1),1)/sqrt(10-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
%plot(squeeze(Drop_ResultCNC_n(1,:,3)),squeeze(nanmean(Drop_ResultCNC_n(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_ResultCNC(11,:,3)),squeeze(nanmean(Drop_ResultCNC(27:end,:,1),1)),squeeze(nanstd(Drop_ResultCNC(27:end,:,1),1)/sqrt(6-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_n=[squeeze(nanmean(Drop_ResultCNC(17:26,:,1),1)),squeeze(nanmean(Drop_ResultCNC(27:end,:,1),1))];
drop_all_expdif=[squeeze(nanmean(Drop_ResultCNC(17:26,:,3),1)),squeeze(nanmean(Drop_ResultCNC(27:end,:,3),1))];

errorbar(squeeze(Drop_ResultCNC(1,:,3)),squeeze(nanmean(Drop_ResultCNC(17:26,:,2),1)),squeeze(nanstd(Drop_ResultCNC(17:26,:,2),1)/sqrt(10-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
errorbar(squeeze(Drop_ResultCNC(11,:,3)),squeeze(nanmean(Drop_ResultCNC(27:end,:,2),1)),squeeze(nanstd(Drop_ResultCNC(27:end,:,2),1)/sqrt(6-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_i=[squeeze(nanmean(Drop_ResultCNC(17:26,:,2),1)),squeeze(nanmean(Drop_ResultCNC(27:end,:,2),1))];
drop_all_expdif=[squeeze(nanmean(Drop_ResultCNC(17:26,:,3),1)),squeeze(nanmean(Drop_ResultCNC(27:end,:,3),1))];

exp_def = fittype('(a*exp(-a*(x)))');
II=~isnan(drop_all_n);
f_1 = fit(drop_all_expdif(II)',drop_all_n(II)',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1, 0]);%,'Lower', [0, -inf, 0, -inf]);
[~,aa_n]=predint(f_1,-5:0.1:7);
coeffvalues(f_1)
plot(-5:0.1:7,aa_n,'Color',[0, 0.5,0],'linewidth',1);hold on;
f_2 = fit(drop_all_expdif(~isnan(drop_all_i))',drop_all_i(~isnan(drop_all_i))',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1,  0]);%,,'Lower', [0, -inf, 0, -inf]);
[~,aa_i]=predint(f_2,-5:0.1:7);
coeffvalues(f_2)
plot(-5:0.1:7,aa_i,'Color',[1, 0.5,0],'linewidth',1)
axis([-6 8 -0.1 0.8])
box off;set(gca,'TickDir','out')
axis square




exp_def = fittype('(a*exp(-a*(x)))');

for ii=1:size(Drop_Result)
    a=isnan(Drop_ResultCNC(ii,:,1))==0 & isnan(squeeze(Drop_ResultCNC(ii,:,2)))==0;
    f1 = fit(squeeze(Drop_ResultCNC(ii,a,3))',squeeze(Drop_ResultCNC(ii,a,1))',exp_def,'Robust','LAR');%,'Robust','LAR');%,'Upper', [10, -8],'Lower', [0.5, -20]);
    f2 = fit(squeeze(Drop_ResultCNC(ii,a,3))',squeeze(Drop_ResultCNC(ii,a,2))',exp_def,'Robust','LAR');%,'Robust','LAR');%,'Upper', [10, -8],'Lower', [0.5, -20]);

    Drop_Result_con_coef(ii,:)=coeffvalues(f1);
    Drop_Result_inact_coef(ii,:)=coeffvalues(f2);

end
Drop_result_con_coef_choice=Drop_Result_con_coef;
Drop_result_inact_coef_choice=Drop_Result_inact_coef;


% Drop_Result_con_coef(Drop_Result_con_coef<-24.5)=nan;
% Drop_Result_inact_coef(Drop_Result_inact_coef<-24.5)=nan;

subplot(243)
a=[nanmean(Drop_Result_con_coef(1:16,1)),nanmean(Drop_Result_inact_coef(1:16,1));...
    nanmean(Drop_Result_con_coef(17:end,1)),nanmean(Drop_Result_inact_coef(17:end,1))];
a_nanstd=[nanstd(Drop_Result_con_coef(1:16,1))/sqrt(16-1),nanstd(Drop_Result_inact_coef(1:16,1))/sqrt(16-1);...
    nanstd(Drop_Result_con_coef(17:end,1))/sqrt(15-1),nanstd(Drop_Result_inact_coef(17:end,1))/sqrt(15-1)];
errorbarplot(a,a_nanstd);
box off;set(gca,'TickDir','out')
axis([0.5 2.5 0 0.25]);
axis square


subplot(244)
plot((ones(16,1)*[1,2])',([Drop_Result_con_coef(1:16,1),Drop_Result_inact_coef(1:16,1)])','k','MarkerSize',15);hold on
plot((ones(15,1)*[3,4])',([Drop_Result_con_coef(17:end,1),Drop_Result_inact_coef(17:end,1)])','k','MarkerSize',15);hold on
plot(ones(1,16),Drop_Result_con_coef(1:16,1),'.','Color',[0, 0.5,0],'MarkerSize',15);hold on
plot(2*ones(1,16),Drop_Result_inact_coef(1:16,1),'.','Color',[1, 0.5,0],'MarkerSize',15);hold on
plot(3*ones(1,15),Drop_Result_con_coef(17:end,1),'.','Color',[0, 0.5,0],'MarkerSize',15);hold on
plot(4*ones(1,15),Drop_Result_inact_coef(17:end,1),'.','Color',[1, 0.5,0],'MarkerSize',15);hold on

[~,p(3)]=ttest(Drop_Result_con_coef,Drop_Result_inact_coef)
[~,p(1)]=ttest(Drop_Result_con_coef(1:16,1),Drop_Result_inact_coef(1:16,1));
[~,p(2)]=ttest(Drop_Result_con_coef(17:end,1),Drop_Result_inact_coef(17:end,1));
axis square
box off;set(gca,'TickDir','out')
axis([0 5 0 0.3]);


subplot(245)
%plot(squeeze(Drop_Result_g(1,:,3)),squeeze(nanmean(Drop_Result_g(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_ResultCNC(1,:,6)),squeeze(nanmean(Drop_ResultCNC(1:10,:,4),1)),squeeze(nanstd(Drop_ResultCNC(1:10,:,4),1)/sqrt(10-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
%plot(squeeze(Drop_ResultCNC_n(1,:,3)),squeeze(nanmean(Drop_ResultCNC_n(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_ResultCNC(11,:,6)),squeeze(nanmean(Drop_ResultCNC(11:16,:,4),1)),squeeze(nanstd(Drop_ResultCNC(11:16,:,4),1)/sqrt(6-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_n=[squeeze(nanmean(Drop_ResultCNC(1:10,:,4),1)),squeeze(nanmean(Drop_ResultCNC(11:16,:,4),1))];
drop_all_expdif=[squeeze(nanmean(Drop_ResultCNC(1:10,:,6),1)),squeeze(nanmean(Drop_ResultCNC(11:16,:,6),1))];

errorbar(squeeze(Drop_ResultCNC(1,:,6)),squeeze(nanmean(Drop_ResultCNC(1:10,:,5),1)),squeeze(nanstd(Drop_ResultCNC(1:10,:,5),1)/sqrt(10-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
errorbar(squeeze(Drop_ResultCNC(11,:,6)),squeeze(nanmean(Drop_ResultCNC(11:16,:,5),1)),squeeze(nanstd(Drop_ResultCNC(11:16,:,5),1)/sqrt(6-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_i=[squeeze(nanmean(Drop_ResultCNC(1:10,:,5),1)),squeeze(nanmean(Drop_ResultCNC(11:16,:,5),1))];
drop_all_expdif=[squeeze(nanmean(Drop_ResultCNC(1:10,:,6),1)),squeeze(nanmean(Drop_ResultCNC(11:16,:,6),1))];

exp_def = fittype('(a*exp(-a*(x)))');
II=~isnan(drop_all_n);
f_1 = fit(drop_all_expdif(II)',drop_all_n(II)',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1, 0]);%,'Lower', [0, -inf, 0, -inf]);
[~,aa_n]=predint(f_1,-5:0.1:7);
coeffvalues(f_1)
plot(-5:0.1:7,aa_n,'Color',[0, 0.5,0],'linewidth',1);hold on;
f_2 = fit(drop_all_expdif(~isnan(drop_all_i))',drop_all_i(~isnan(drop_all_i))',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1,  0]);%,,'Lower', [0, -inf, 0, -inf]);
[~,aa_i]=predint(f_2,-5:0.1:7);
coeffvalues(f_2)
plot(-5:0.1:7,aa_i,'Color',[1, 0.5,0],'linewidth',1)
axis([-6 8 -0.1 1])
box off;set(gca,'TickDir','out')
axis square


subplot(246)
%plot(squeeze(Drop_ResultCNC_g(1,:,3)),squeeze(nanmean(Drop_ResultCNC_g(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_ResultCNC(1,:,6)),squeeze(nanmean(Drop_ResultCNC(17:26,:,4),1)),squeeze(nanstd(Drop_ResultCNC(17:26,:,4),1)/sqrt(10-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
%plot(squeeze(Drop_ResultCNC_n(1,:,3)),squeeze(nanmean(Drop_ResultCNC_n(:,:,1),1)),'Color',[0, 0.5,0],'Marker', '.','MarkerSize',21);hold on;%%
errorbar(squeeze(Drop_ResultCNC(11,:,6)),squeeze(nanmean(Drop_ResultCNC(27:end,:,4),1)),squeeze(nanstd(Drop_ResultCNC(27:end,:,4),1)/sqrt(6-1)),'Color',[0, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_n=[squeeze(nanmean(Drop_ResultCNC(17:26,:,4),1)),squeeze(nanmean(Drop_ResultCNC(27:end,:,4),1))];
drop_all_expdif=[squeeze(nanmean(Drop_ResultCNC(17:26,:,6),1)),squeeze(nanmean(Drop_ResultCNC(27:end,:,6),1))];

errorbar(squeeze(Drop_ResultCNC(1,:,6)),squeeze(nanmean(Drop_ResultCNC(17:26,:,5),1)),squeeze(nanstd(Drop_ResultCNC(17:26,:,5),1)/sqrt(10-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
errorbar(squeeze(Drop_ResultCNC(11,:,6)),squeeze(nanmean(Drop_ResultCNC(27:end,:,5),1)),squeeze(nanstd(Drop_ResultCNC(27:end,:,5),1)/sqrt(6-1)),'Color',[1, 0.5,0],'LineStyle','none','Marker', '.','MarkerSize',21);hold on;
drop_all_i=[squeeze(nanmean(Drop_ResultCNC(17:26,:,5),1)),squeeze(nanmean(Drop_ResultCNC(27:end,:,5),1))];
drop_all_expdif=[squeeze(nanmean(Drop_ResultCNC(17:26,:,6),1)),squeeze(nanmean(Drop_ResultCNC(27:end,:,6),1))];

exp_def = fittype('(a*exp(-a*(x)))');
II=~isnan(drop_all_n);
f_1 = fit(drop_all_expdif(II)',drop_all_n(II)',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1, 0]);%,'Lower', [0, -inf, 0, -inf]);
[~,aa_n]=predint(f_1,-5:0.1:7);
coeffvalues(f_1)
plot(-5:0.1:7,aa_n,'Color',[0, 0.5,0],'linewidth',1);hold on;
f_2 = fit(drop_all_expdif(~isnan(drop_all_i))',drop_all_i(~isnan(drop_all_i))',exp_def,'Robust','LAR');%,'Upper', [Inf, -0.1,  0]);%,,'Lower', [0, -inf, 0, -inf]);
[~,aa_i]=predint(f_2,-5:0.1:7);
coeffvalues(f_2)
plot(-5:0.1:7,aa_i,'Color',[1, 0.5,0],'linewidth',1)
axis([-6 8 -0.1 0.8])
box off;set(gca,'TickDir','out')
axis square




exp_def = fittype('(a*exp(-a*(x)))');

for ii=1:size(Drop_Result)

    
    a=isnan(Drop_ResultCNC(ii,:,4))==0 & isnan(squeeze(Drop_ResultCNC(ii,:,5)))==0;
    f1 = fit(squeeze(Drop_ResultCNC(ii,a,6))',squeeze(Drop_ResultCNC(ii,a,4))',exp_def,'Robust','LAR');%,'Upper', [10, -8],'Lower', [0.5, -20]);
    f2 = fit(squeeze(Drop_ResultCNC(ii,a,6))',squeeze(Drop_ResultCNC(ii,a,5))',exp_def,'Robust','LAR');%,'Upper', [10, -8],'Lower', [0.5, -20]);

    Drop_Result_con_coef(ii,:)=coeffvalues(f1);
    Drop_Result_inact_coef(ii,:)=coeffvalues(f2);

end
Drop_result_con_coef_nochoice=Drop_Result_con_coef;
Drop_result_inact_coef_nochoice=Drop_Result_inact_coef;


[~,p]=ttest(Drop_result_con_coef_choice(1:16)-Drop_result_con_coef_nochoice(1:16))

nanmean(Drop_result_con_coef_choice(1:16)-Drop_result_con_coef_nochoice(1:16))
ttest(Drop_result_con_coef_choice(1:16)-Drop_result_con_coef_nochoice(1:16))
nanmean(Drop_result_con_coef_choice(17:end)-Drop_result_con_coef_nochoice(17:end))
[~,p]=ttest(Drop_result_con_coef_choice(17:end)-Drop_result_con_coef_nochoice(17:end))

nanmean(Drop_result_con_coef_choice(1:16)-Drop_result_inact_coef_choice(1:16))
ttest(Drop_result_con_coef_choice(1:16)-Drop_result_inact_coef_choice(1:16))

nanmean(Drop_result_con_coef_nochoice-Drop_result_inact_coef_nochoice)
[~,p]=ttest(Drop_result_con_coef_nochoice-Drop_result_inact_coef_nochoice)
nanmean(Drop_result_con_coef_nochoice(1:16)-Drop_result_inact_coef_nochoice(1:16))
[~,p]=ttest(Drop_result_con_coef_nochoice(1:16)-Drop_result_inact_coef_nochoice(1:16))
nanmean(Drop_result_con_coef_nochoice(17:end)-Drop_result_inact_coef_nochoice(17:end))
[~,p]=ttest(Drop_result_con_coef_nochoice(17:end)-Drop_result_inact_coef_nochoice(17:end))

nanmean(Drop_result_con_coef_choice-Drop_result_inact_coef_choice)
[~,p]=ttest(Drop_result_con_coef_choice-Drop_result_inact_coef_choice)
nanmean(Drop_result_con_coef_choice(1:16)-Drop_result_inact_coef_choice(1:16))
[~,p]=ttest(Drop_result_con_coef_choice(1:16)-Drop_result_inact_coef_choice(1:16))
nanmean(Drop_result_con_coef_choice(17:end)-Drop_result_inact_coef_choice(17:end))
[~,p]=ttest(Drop_result_con_coef_choice(17:end)-Drop_result_inact_coef_choice(17:end))

% Drop_Result_con_coef(Drop_Result_con_coef<-24.5)=nan;
% Drop_Result_inact_coef(Drop_Result_inact_coef<-24.5)=nan;

subplot(247)
a=[nanmean(Drop_Result_con_coef(1:16,1)),nanmean(Drop_Result_inact_coef(1:16,1));...
    nanmean(Drop_Result_con_coef(17:end,1)),nanmean(Drop_Result_inact_coef(17:end,1))];
a_nanstd=[nanstd(Drop_Result_con_coef(1:16,1))/sqrt(16-1),nanstd(Drop_Result_inact_coef(1:16,1))/sqrt(16-1);...
    nanstd(Drop_Result_con_coef(17:end,1))/sqrt(15-1),nanstd(Drop_Result_inact_coef(17:end,1))/sqrt(15-1)];
errorbarplot(a,a_nanstd);
box off;set(gca,'TickDir','out')
axis([0.5 2.5 0 0.25]);
[~,p(3)]=ttest(Drop_Result_con_coef,Drop_Result_inact_coef)
[~,p(1)]=ttest(Drop_Result_con_coef(1:16,1),Drop_Result_inact_coef(1:16,1));
[~,p(2)]=ttest(Drop_Result_con_coef(17:end,1),Drop_Result_inact_coef(17:end,1));
axis square


subplot(248)
plot((ones(16,1)*[1,2])',([Drop_Result_con_coef(1:16,1),Drop_Result_inact_coef(1:16,1)])','k','MarkerSize',15);hold on
plot((ones(15,1)*[3,4])',([Drop_Result_con_coef(17:end,1),Drop_Result_inact_coef(17:end,1)])','k','MarkerSize',15);hold on

plot(ones(1,16),Drop_Result_con_coef(1:16,1),'.','Color',[0, 0.5,0],'MarkerSize',15);hold on
plot(2*ones(1,16),Drop_Result_inact_coef(1:16,1),'.','Color',[1, 0.5,0],'MarkerSize',15);hold on
plot(3*ones(1,15),Drop_Result_con_coef(17:end,1),'.','Color',[0, 0.5,0],'MarkerSize',15);hold on
plot(4*ones(1,15),Drop_Result_inact_coef(17:end,1),'.','Color',[1, 0.5,0],'MarkerSize',15);hold on
box off;set(gca,'TickDir','out')
axis([0 5 0 0.3]);
[~,p(3)]=ttest(Drop_Result_con_coef,Drop_Result_inact_coef)
[~,p(1)]=ttest(Drop_Result_con_coef(1:16,1),Drop_Result_inact_coef(1:16,1));
[~,p(2)]=ttest(Drop_Result_con_coef(17:end,1),Drop_Result_inact_coef(17:end,1));
axis square
%%
I_A_1=[ones(1,16),zeros(1,15)]'==1 & sum(FixControlH1d,2)~=0;
I_I_1=[zeros(1,16),ones(1,15)]'==1 &  sum(FixControlH1d,2)~=0;
I_A_2=[ones(1,16),zeros(1,15)]'==1 & sum(FixInactH1d,2)~=0;
I_I_2=[zeros(1,16),ones(1,15)]'==1 &  sum(FixInactH1d,2)~=0;


colormap hot
clim=[0 0.0005];
figure(1)
subplot(231)
%a=imgaussfilt(squeeze(nanmean(FixControlH,1)),0.1);
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(FixControlH(:,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(232)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(SacControlH(:,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(233)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(SucControlH(:,:,:),1))+squeeze(nanmean(FaiControlH(:,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(234)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(FixInactH(:,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(235)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(SacInactH(:,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(236)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(SucInactH(:,:,:),1))+squeeze(nanmean(FaiInactH(:,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);
clim=[0 0.001];



colormap hot
clim=[0 0.0005];
figure(2)
subplot(231)
%a=imgaussfilt(squeeze(nanmean(FixControlH,1)),0.1);
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(FixControlH(I_I_1,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(232)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(SacControlH(I_I_1,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(233)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(SucControlH(I_I_1,:,:),1))+squeeze(nanmean(FaiControlH(17:end,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(234)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(FixInactH(I_I_2,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(235)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(SacInactH(I_I_2,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);

subplot(236)
[a,xbins0,ybins0]=interpretHist(squeeze(nanmean(SucInactH(I_I_2,:,:),1))+squeeze(nanmean(FaiInactH(17:end,:,:),1)),xbins,ybins);
imagesc(xbins,ybins,a,clim);  box off;set(gca,'TickDir','out');axis square
axis([-25 25 -25 25]);
clim=[0 0.001];
colormap hot
%         subplot(244)
%         imagesc(xbins,ybins,squeeze(nanmean(FaiControlH,1)),clim);  box off;set(gca,'TickDir','out');axis square
%         subplot(248)
%         imagesc(xbins,ybins,squeeze(nanmean(FaiInactH,1)),clim); box off;set(gca,'TickDir','out');axis square

figure(3)

subplot(231)
plotstd(xbins,(SacControlH1d(I_A_1,:)),'r');hold on;
plotstd(xbins,(SacInactH1d(I_A_2,:)),'b');hold on;
[~,p]=kstest2(nanmean(SacControlH1d(I_A_1,:),1),nanmean(SacInactH1d(I_A_2,:),1))
axis([5 15 0 0.6]); title(['p=',num2str(p,3)]);
box off;set(gca,'TickDir','out');axis square

subplot(234)
plotstd(xbins,SacControlH1d(I_I_1,:),'r');hold on;
plotstd(xbins,SacInactH1d(I_I_2,:),'b');hold on;
[~,p]=kstest2(nanmean(SacControlH1d(I_I_1,:),1),nanmean(SacInactH1d(I_I_2,:),1))
axis([5 15 0 0.6]); title(['p=',num2str(p,3)]);
box off;set(gca,'TickDir','out');axis square

subplot(232)
plotstd(xbins,(SucControlH1d(I_A_1,:)),'r');hold on;
plotstd(xbins,(SucInactH1d(I_A_2,:)),'b');hold on;
[~,p]=kstest2(nanmean(SucControlH1d(I_A_1,:),1),nanmean(SucInactH1d(I_A_2,:),1))
axis([5 15 0 0.6]); title(['p=',num2str(p,3)]);
box off;set(gca,'TickDir','out');axis square

subplot(235)
plotstd(xbins,(SucControlH1d(I_I_1,:)),'r');hold on;
plotstd(xbins,(SucInactH1d(I_I_2,:)),'b');hold on;
[~,p]=kstest2(nanmean(SucControlH1d(1:16,:),1),nanmean(SucInactH1d(1:16,:),1))
axis([5 15 0 1]); title(['p=',num2str(p,3)]);
box off;set(gca,'TickDir','out');axis square
% subplot(233)
% per
% bar([Break_sac

% % subplot(233)
% % a= sum(Break_sac(I_A_1,:)==1,2)./(sum(Break_sac(I_A_1,:)==1,2) + sum(Break_sac(I_A_1,:)==-1,2));
% % b= sum(Break_sac(I_A_1,:)==2,2)./(sum(Break_sac(I_A_1,:)==2,2) + sum(Break_sac(I_A_1,:)==-2,2));
% % bar([nanmean(a),nanmean(b)],0.4);hold on;
% % errorbar([nanmean(a),nanmean(b)],[nanstd(a),nanstd(b)]./sqrt([sum(~isnan(a))-1, sum(~isnan(b))-1]),'k.');hold on;
% % axis([0 3 0 1])
% % [~,p1]=ttest(a,b);
% % title(['p=',num2str(p1)]);
% % box off;set(gca,'TickDir','out');axis square
% %
% % subplot(236)
% % a= sum(Break_sac(I_I_1,:)==1,2)./(sum(Break_sac(I_I_1,:)==1,2) + sum(Break_sac(I_I_1,:)==-1,2));
% % b= sum(Break_sac(I_I_1,:)==2,2)./(sum(Break_sac(I_I_1,:)==2,2) + sum(Break_sac(I_I_1,:)==-2,2));
% % bar([nanmean(a),nanmean(b)],0.4);hold on;
% % errorbar([nanmean(a),nanmean(b)],[nanstd(a),nanstd(b)]./sqrt([sum(~isnan(a))-1, sum(~isnan(b))-1]),'k.');hold on;
% % axis([0 3 0 1])
% % [~,p1]=ttest(a,b);
% % title(['p=',num2str(p1)]);
% % box off;set(gca,'TickDir','out');axis square

%%
figure(4)
subplot(221)
errorbarplot([nanmean(FixateControlDelta(1:16)),nanmean(SacControlDelta(1:16)),nanmean(SucControlDelta(1:16));...
    nanmean(FixateInactDelta(1:16)),nanmean(SacInactDelta(1:16)),nanmean(SucInactDelta(1:16))]',...
    [nanstderror(FixateControlDelta(1:16)),nanstderror(SacControlDelta(1:16)),nanstderror(SucControlDelta(1:16));...
    nanstderror(FixateInactDelta(1:16)),nanstderror(SacInactDelta(1:16)),nanstderror(SucInactDelta(1:16))]',0.8 );
box off;set(gca,'TickDir','out');axis square
[h,p1]=ttest(FixateControlDelta(1:16),FixateInactDelta(1:16))
[h,p2]=ttest(SacControlDelta(1:16),SacInactDelta(1:16))
[h,p3]=ttest(SucControlDelta(1:16),SucInactDelta(1:16))
axis([0 4 0 1.5]);

subplot(222)
errorbarplot([nanmean(FixateControlDelta(17:end)),nanmean(SacControlDelta(17:end)),nanmean(SucControlDelta(17:end));...
    nanmean(FixateInactDelta(17:end)),nanmean(SacInactDelta(17:end)),nanmean(SucInactDelta(17:end))]',...
    [nanstderror(FixateControlDelta(17:end)),nanstderror(SacControlDelta(17:end)),nanstderror(SucControlDelta(17:end));...
    nanstderror(FixateInactDelta(17:end)),nanstderror(SacInactDelta(17:end)),nanstderror(SucInactDelta(17:end))]',0.8 );

[h,p1]=ttest(FixateControlDelta(17:end),FixateInactDelta(17:end))
[h,p2]=ttest(SacControlDelta(17:end),SacInactDelta(17:end))
[h,p3]=ttest(SucControlDelta(17:end),SucInactDelta(17:end))
box off;set(gca,'TickDir','out');axis square
axis([0 4 0 1.5]);


subplot(223)
errorbarplot([nanmean(V_S(1:16,:),1);nanmean(V_S(17:end,:),1)],[nanstderror(V_S(1:16,:),1);nanstderror(V_S(17:end,:),1)],0.8);
box off;set(gca,'TickDir','out');axis square
[h,p1]=ttest(V_S(1:16,1),V_S(1:16,2))
[h,p2]=ttest(V_S(17:end,1),V_S(17:end,2))
axis([0 3 0 450]);


subplot(224)
errorbarplot([nanmean(ResultBreak(1:16,:),1);nanmean(ResultBreak(17:end,:),1)],[nanstderror(ResultBreak(1:16,:),1);nanstderror(ResultBreak(17:end,:),1)],0.8);
box off;set(gca,'TickDir','out');axis square
[h,p1]=ttest(ResultBreak(1:16,1),ResultBreak(1:16,2))
[h,p2]=ttest(ResultBreak(17:end,1),ResultBreak(17:end,2))




% errorbarplot([nanmean(FixateControlDelta(1:16)),nanmean(FixateControlDelta(17:end)); nanmean(FixateInactDelta(1:16)),nanmean(FixateInactDelta(17:end))]',...
%   [nanstderror(FixateControlDelta(1:16)),nanstderror(FixateControlDelta(17:end)); nanstderror(FixateInactDelta(1:16)),nanstderror(FixateInactDelta(17:end))]  );
% subplot(222)
% errorbarplot([nanmean(SacControlDelta(1:16)),nanmean(SacControlDelta(17:end)); nanmean(SacInactDelta(1:16)),nanmean(SacInactDelta(17:end))]',...
%   [nanstderror(SacControlDelta(1:16)),nanstderror(SacControlDelta(17:end)); nanstderror(SacInactDelta(1:16)),nanstderror(SacInactDelta(17:end))]  );
% subplot(223)
% errorbarplot([nanmean(SucControlDelta(1:16)),nanmean(SucControlDelta(17:end)); nanmean(SucInactDelta(1:16)),nanmean(SucInactDelta(17:end))]',...
%   [nanstderror(SucControlDelta(1:16)),nanstderror(SucControlDelta(17:end)); nanstderror(SucInactDelta(1:16)),nanstderror(SucInactDelta(17:end))]  );

%%
figure(5)
subplot(2,2,1)
errorbarplot(squeeze(nanmean(Sac_setion_noChoice(1:16,:,:),1)),squeeze(nanstd(Sac_setion_noChoice(1:16,:,:),1))/sqrt(15));
box off;set(gca,'TickDir','out')
axis([0.5 4.5 0 220]);

subplot(2,2,2)
errorbarplot(squeeze(nanmean(Sac_setion_Choice(1:16,:,:),1)),squeeze(nanstd(Sac_setion_Choice(1:16,:,:),1))/sqrt(15));
box off;set(gca,'TickDir','out')
axis([0.5 4.5 0 220]);

subplot(2,2,3)
errorbarplot(squeeze(nanmean(Sac_setion_noChoice(17:end,:,:),1)),squeeze(nanstd(Sac_setion_noChoice(17:end,:,:),1))/sqrt(14));
box off;set(gca,'TickDir','out')
axis([0.5 4.5 0 220]);

subplot(2,2,4)
errorbarplot(squeeze(nanmean(Sac_setion_Choice(17:end,:,:),1)),squeeze(nanstd(Sac_setion_Choice(17:end,:,:),1))/sqrt(14));
box off;set(gca,'TickDir','out')
axis([0.5 4.5 0 220]);


a=reshape(Sac_setion_noChoice(1:16,:,:),[64 2]);
anova2(a,16)


a=reshape(Sac_setion_Choice(1:16,:,:),[64 2]);
anova2(a,16)


a=reshape(Sac_setion_noChoice(17:end,:,:),[60 2]);
anova2(a,15)

a=reshape(Sac_setion_Choice(17:end,:,:),[60 2]);
anova2(a,15)


% subplot(3,2,5)
% plot(squeeze(nanmean(Sac_setion_noChoice(1:16,:,1),2)),squeeze(nanmean(Sac_setion_noChoice(1:16,:,2),2)),'k.','MarkerSize',10);hold on;
% plot(squeeze(nanmean(Sac_setion_noChoice(17:end,:,1),2)),squeeze(nanmean(Sac_setion_noChoice(17:end,:,2),2)),'b.','MarkerSize',10);hold on;
%
% plot(100:1:220,100:1:220,'k');
% axis([140 220 140 220]);
% axis square
% [h,p]=ttest(squeeze(nanmean(Sac_setion_noChoice(:,:,1),2)),squeeze(nanmean(Sac_setion_noChoice(:,:,2),2)));
% title(['p=',num2str(p,3)]);
% box off;set(gca,'TickDir','out');
%
% subplot(3,2,6)
% plot(squeeze(nanmean(Sac_setion_Choice(1:16,:,1),2)),squeeze(nanmean(Sac_setion_Choice(1:16,:,2),2)),'k.','MarkerSize',10);hold on;
% plot(squeeze(nanmean(Sac_setion_Choice(17:end,:,1),2)),squeeze(nanmean(Sac_setion_Choice(17:end,:,2),2)),'b.','MarkerSize',10);hold on;
%
% plot(100:1:220,100:1:220,'k');
% axis([140 220 140 220]);
% axis square
% [h,p]=ttest(squeeze(nanmean(Sac_setion_Choice(:,:,1),2)),squeeze(nanmean(Sac_setion_Choice(:,:,2),2)));
% title(['p=',num2str(p,3)]);
% box off;set(gca,'TickDir','out');

%%
figure(6)
V_S(V_S==0)=nan;
V_S=reshape(V_S,[31,2]);
T_all(T_all==0)=nan;
T_all=reshape(T_all,[31,2]);
I_A=1:15;
subplot(221)
[h,p]=ttest(V_S(I_A,1),V_S(I_A,2))
bar(nanmean(V_S(I_A,:),1),0.4);hold on;
errorbar(1:2,nanmean(V_S(I_A,:),1),nanstd(V_S(I_A,:),1)/sqrt(14),'.k');
axis([0 3 200 400]);
box off;set(gca,'TickDir','out');
title(['p=',num2str(p)]);

subplot(222)
[h,p]=ttest(T_all(I_A,1),T_all(I_A,2))
bar(nanmean(T_all(I_A,:),1),0.4);hold on;
errorbar(1:2,nanmean(T_all(I_A,:),1),nanstd(T_all(I_A,:),1)/sqrt(14),'.k');
axis([0 3 100 300]);
box off;set(gca,'TickDir','out');
title(['p=',num2str(p)]);


subplot(223)
[h,p]=ttest(Dif_fix_x(I_A,1),Dif_fix_x(I_A,2))
bar(nanmean(Dif_fix_x(I_A,:),1),0.4);hold on;
errorbar(1:2,nanmean(Dif_fix_x(I_A,:),1),nanstd(Dif_fix_x(I_A,:),1)/sqrt(14),'.k');
box off;set(gca,'TickDir','out');
axis([0 3 0 0.5]);
title(['p=',num2str(p)]);

subplot(224)
[h,p]=ttest(Dif_fix_y(I_A,1),Dif_fix_y(I_A,2))
bar(nanmean(Dif_fix_y(I_A,:),1),0.4);hold on;
errorbar(1:2,nanmean(Dif_fix_y(I_A,:),1),nanstd(Dif_fix_y(I_A,:),1)/sqrt(14),'.k');
box off;set(gca,'TickDir','out');
axis([0 3 0 0.5]);
title(['p=',num2str(p)]);

%%
figure(3)
I_I=16:29;
subplot(221)
[h,p]=ttest(V_S(I_I,1),V_S(I_I,2))
bar(nanmean(V_S(I_I,:),1),0.4);hold on;
errorbar(1:2,nanmean(V_S(I_I,:),1),nanstd(V_S(I_I,:),1)/sqrt(13),'.k');
axis([0 3 200 400]);
box off;set(gca,'TickDir','out');
title(['p=',num2str(p)]);

subplot(222)
[h,p]=ttest(T_all(I_I,1),T_all(I_I,2))
bar(nanmean(T_all(I_I,:),1),0.4);hold on;
errorbar(1:2,nanmean(T_all(I_I,:),1),nanstd(T_all(I_I,:),1)/sqrt(13),'.k');
axis([0 3 100 300]);
box off;set(gca,'TickDir','out');
title(['p=',num2str(p)]);


subplot(223)
[h,p]=ttest(Dif_fix_x(I_I,1),Dif_fix_x(I_I,2))
bar(nanmean(Dif_fix_x(I_I,:),1),0.4);hold on;
errorbar(1:2,nanmean(Dif_fix_x(I_I,:),1),nanstd(Dif_fix_x(I_I,:),1)/sqrt(13),'.k');
box off;set(gca,'TickDir','out');
axis([0 3 0 0.4]);
title(['p=',num2str(p)]);

subplot(224)
[h,p]=ttest(Dif_fix_y(I_I,1),Dif_fix_y(I_I,2))
bar(nanmean(Dif_fix_y(I_I,:),1),0.4);hold on;
errorbar(1:2,nanmean(Dif_fix_y(I_I,:),1),nanstd(Dif_fix_y(I_I,:),1)/sqrt(13),'.k');
box off;set(gca,'TickDir','out');
axis([0 3 0 0.4]);
title(['p=',num2str(p)]);

h10=figure(2);
print( h10, '-djpeg', ['D:\Projects\SEFCooling\figures\EyeSum_monkeyA']);
print( h10, '-depsc', ['D:\Projects\SEFCooling\figures\EyeSum_monkeyA']);
h10=figure(3);
print( h10, '-djpeg', ['D:\Projects\SEFCooling\figures\EyeSum_monkeyI']);
print( h10, '-depsc', ['D:\Projects\SEFCooling\figures\EyeSum_monkeyI']);
axis([0 3 200 400]);

% for per=1:1000
%     
%     z_x_0=z_x(  randperm(length(z_x)));
%     
%     z_y_0=z_y(  randperm(length(z_x)));
%     
%     
%     
%     dif_fix_x_perm(per)=nanmean(z_x_0(1:length(fix1_x)))-nanmean(z_x_0((length(fix1_x)+1):end));
%     
%     dif_fix_y_perm(per)=nanmean(z_y_0(1:length(fix1_y)))-nanmean(z_y_0((length(fix1_y)+1):end));
%     
% end



%
