function TimeFrequencyPlot_normal(Dirfigure,name,file,eventfile)
close all;
load(file,'Map','ChannelNumber');
load(eventfile);


t_MP=Map.t;
f_MP=Map.f;
clim=[-1 1];
figure(1)
for ch=1:ChannelNumber
    subplot(4,2,1+(ch-1)*2)
    eval(['a=Map.TargMean',num2str(ch),';']);
    tfplot(t_MP,f_MP,a,clim)
    %  imagesc(squeeze(nanmean(Map.Targ1,1)))
    axis([-100 400 0 200]);
    % figure(gcf);set(gca,'YDir','norm');
    subplot(4,2,2+(ch-1)*2)
    eval(['a=Map.ResultMean',num2str(ch),';']);
    tfplot(t_MP,f_MP,a,clim)
    %  imagesc(squeeze(nanmean(Map.Targ1,1)))
    axis([-300 400 0 200]);
end
h=figure(1)
print( h, '-djpeg', [Dirfigure,name,'MP']);

%%

clim=[-2 2];
Targ_t=downsample(t_MP(t_MP>-100 & t_MP<400),10);
Result_t=downsample(t_MP(t_MP>-400 & t_MP<400),10);
Targ_f=downsample(f_MP,2);


for ch=1:ChannelNumber
    figure(ch+1)
    clear a b
    eval(['a=Map.Targ',num2str(ch),';']);
    eval(['b=Map.Result',num2str(ch),';']);
    
    for dir=1:4
        subplot(3,4,dir)
        clear I_dir a1_1
        
        I_dir=TrialInfo{:,15}==dir ;
        a1_1=squeeze(nanmean(a(I_dir,:,:),1));
        tfplot(Targ_t,Targ_f,a1_1,clim)
        axis([-100 400 0 200]);
    end
    
    
    u=0;
    for r=[1 3 5 9]
        u=u+1;
        subplot(3,4,4+u)
        clear I_dir a1_1
        I_dir=TrialInfo{:,17}==r ;
        a1_1=squeeze(nanmean(b(I_dir,:,:),1));
        tfplot(Result_t,Targ_f,a1_1,clim)
        %  imagesc(squeeze(nanmean(Map.Targ1,1)))
        axis([-300 400 0 200]);
    end
    h=figure(ch+1);
    print( h, '-djpeg', [Dirfigure,name,'Ch',num2str(ch),'D_R']);
    
    
    % figure(gcf);set(gca,'YDir','norm');
    u=0;
    for r=[1 3 5 9]
        u=u+1;
        subplot(3,4,4+u)
        clear I_dir a1_1
        I_dir=TrialInfo{:,17}==r ;
        a1_1=squeeze(nanmean(b(I_dir,:,:),1));
        tfplot(Result_t,Targ_f,a1_1,clim)
        %  imagesc(squeeze(nanmean(Map.Targ1,1)))
        axis([-300 400 0 200]);
    end
    h=figure(ch+1);
    print( h, '-djpeg', [Dirfigure,name,'Ch',num2str(ch),'D_R']);

end


