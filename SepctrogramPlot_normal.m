load(file,'Map');

clim=[-1 1];
figure(1)
for ch=1:ch_num
subplot(3,2,1+(ch-1)*2)
eval(['a=Map.TargMean',num2str(ch),';']);
tfplot(t_MP,f_MP,a,clim)
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-100 400 0 200]);
% figure(gcf);set(gca,'YDir','norm');
subplot(3,2,2+(ch-1)*2)
eval(['a=Map.ResultMean',num2str(ch),';']);
tfplot(t_MP,f_MP,a,clim)
%  imagesc(squeeze(nanmean(Map.Targ1,1)))
axis([-300 400 0 200]);
end