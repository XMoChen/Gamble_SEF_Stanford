                subplot(221)
                
                if d==2 &  temp==1 & sum(I0)>5
                    
                    V0_Sac_1=[V0_Sac_1,V0_Sac];
                    T_Sac_1=[T_Sac_1,T_Sac];
                    
                    XX=sqrt(X0_Sac(I0,:).^2+Y0_Sac(I0,:).^2);
                    plot(-500:900,XX(1:5,:),'r');hold on
                    
                elseif d==4 &  sum(I0)>5
                    
                    V0_Sac_2=[V0_Sac_2,V0_Sac];
                    T_Sac_2=[T_Sac_2,T_Sac(I0)];
                    
                    XX=sqrt(X0_Sac(I0,:).^2+Y0_Sac(I0,:).^2);
                    plot(-500:900,XX(1:5,:),'b');hold on
                end
                
                axis([-300 150 -5 25]);
                plot(-300:0,2*ones(1,301),'k--');hold on;
                plot(-300:0,-2*ones(1,301),'k--');hold on;
                plot(0:100,(10+3)*ones(1,101),'k--');hold on;
                plot(0:100,(10-3)*ones(1,101),'k--');hold on;
                plot(0:100,(-10+3)*ones(1,101),'k--');hold on;
                plot(0:100,(-10-3)*ones(1,101),'k--');
                
                box off;set(gca,'TickDir','out');
                
                axis square
                
                
                
                subplot(222)
                
                if d==1 &  temp==1 & sum(I1)>5
                    YY=sqrt(X0_Sac(I1,:).^2+Y0_Sac(I1,:).^2);
                    plot(-500:900,YY(1:5,:),'r');hold on
                    
                elseif d==4 & sum(I1)>5
                    YY=sqrt(X0_Sac(I1,:).^2+Y0_Sac(I1,:).^2);
                    plot(-500:900,YY(1:5,:),'b');hold on
                    
                end
                
                axis([-300 150 -5 25]);
                plot(-100:0,2*ones(1,101),'k--');hold on;
                plot(-100:0,-2*ones(1,101),'k--');hold on;
                plot(0:100,(10+3)*ones(1,101),'k--');hold on;
                plot(0:100,(10-3)*ones(1,101),'k--');hold on;
                plot(0:100,(-10+3)*ones(1,101),'k--');hold on;
                plot(0:100,(-10-3)*ones(1,101),'k--');
                box off;set(gca,'TickDir','out');
                
                axis square

                
               subplot(223)
                
                if  d==2 &  temp==1 & sum(I0)>5
            
                    
                    XX=sqrt(X_reward(I0,:).^2+Y_reward(I0,:).^2);
                    plot(-300:200,XX(1:6,:),'r');hold on
                    
                elseif d==2 &  sum(I0)>5

                    XX=sqrt(X_reward(I0,:).^2+Y_reward(I0,:).^2);
                    plot(-300:200,XX(1:5,:),'b');hold on
                end
                
                axis([-300 150 -5 25]);
%                 plot(-100:0,2*ones(1,101),'k--');hold on;
%                 plot(-100:0,-2*ones(1,101),'k--');hold on;
                plot(-300:100,(10+3)*ones(1,401),'k--');hold on;
                plot(-300:100,(10-3)*ones(1,401),'k--');hold on;
%                 plot(0:100,(-10+3)*ones(1,101),'k--');hold on;
%                 plot(0:100,(-10-3)*ones(1,101),'k--');
                
                box off;set(gca,'TickDir','out');axis square
                
                
                
                subplot(224)
                
                if d==2 &  temp==1 & sum(I1)>5
                    YY=sqrt(X_reward(I1,:).^2+Y_reward(I1,:).^2);
                    plot(-300:200,YY(1:5,:),'r');hold on
                    
                elseif d==4 &  sum(I1)>5
                    YY=sqrt(X_reward(I1,:).^2+Y_reward(I1,:).^2);
                    plot(-300:200,YY(1:5,:),'b');hold on
                    
                end
                plot(-300:100,(10+3)*ones(1,401),'k--');hold on;
                plot(-300:100,(10-3)*ones(1,401),'k--');hold on;                
                 axis([-300 150 -5 25]);
%                 plot(-100:0,2*ones(1,101),'k--');hold on;
%                 plot(-100:0,-2*ones(1,101),'k--');hold on;
%                 plot(0:100,(10+3)*ones(1,101),'k--');hold on;
%                 plot(0:100,(10-3)*ones(1,101),'k--');hold on;
%                 plot(0:100,(-10+3)*ones(1,101),'k--');hold on;
%                 plot(0:100,(-10-3)*ones(1,101),'k--');
                box off;set(gca,'TickDir','out');
                
                axis square
                
