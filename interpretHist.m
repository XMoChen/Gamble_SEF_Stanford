function [im,Xq,Yq]=interpretHist(a,xtarget,ytarget)
[Xq1,Yq1] = meshgrid(xtarget,ytarget);
    Xq = (min(xtarget)):0.2:(max(xtarget));
    Yq = (min(ytarget)):0.2:(max(ytarget));
    [Xq2,Yq2] = meshgrid(Xq ,Yq);
    im= interp2(Xq1,Yq1,a,Xq2,Yq2);