function [AIC_sum,b1,b2,b3,b4]=modelfitting(Move50,PreferedDirection,table,PsiV)


%%%%% chosen target is in the receptive fild
[xCDir,xNCDir,CValue,NCValue,CPsiValue,NCPsiValue]=ChNChDirection2(table,PsiV);
IC_Pre=xCDir==PreferedDirection;
Vin=CPsiValue(IC_Pre);
Vout=NCPsiValue(IC_Pre);
x(:,1)=Vin;
x(:,2)=Vout;
x(:,3)=table(IC_Pre,23)~=0;
Rmax=max(Move50);
R=Move50(IC_Pre);
iterationno=100;

for j = 1:iterationno
%%fraction model
b0 = [Rmax*rand(1) Rmax*rand(1) Rmax*rand(1) rand(1)];b0=double(b0);
ub = [Rmax  Rmax  Rmax 1];ub=double(ub)';
lb = [0   0 0  0];lb=double(lb)';
options = optimset('lsqnonlin');

myfun1 = @(b) (R-(b(1)+b(2).*(x(:,1))./((x(:,1))+b(3).*(x(:,2)))));
optnew=optimset(options,'MaxFunEvals',10000, 'MaxIter', 3000);
[b1,rss1] = lsqnonlin(myfun1,b0,lb,ub,optnew);

%%difference contribution model
b0 = [Rmax*rand(1) Rmax*rand(1) Rmax*rand(1) rand(1)];b0=double(b0);
ub = [Rmax  Rmax  Rmax 1];ub=double(ub)';
lb = [0   -Rmax  -Rmax 0];lb=double(lb)';

myfun2 = @(b) (R-(b(1)+b(2).*(x(:,1))+b(3).*(x(:,2))));
optnew=optimset(options,'MaxFunEvals',10000, 'MaxIter', 3000);
[b2,rss2] = lsqnonlin(myfun2,b0,lb,ub,optnew);

%%simple divisive normalization model
b0 = [rand(1) rand(1) rand(1)];b0=double(b0);
ub = [10  10  1];ub=double(ub)';
lb = [0    0  0];lb=double(lb)';

myfun3 = @(b) (R-Rmax*((x(:,1))./(b(2)+(x(:,1))+b(1).*(x(:,2)))));
optnew=optimset(options,'MaxFunEvals',10000, 'MaxIter', 3000);
[b3,rss3] = lsqnonlin(myfun3,b0,lb,ub,optnew);


%%divisive normalization model
b0 = [rand(1) rand(1) rand(1)  rand(1) ];b0=double(b0);
ub = [10 10 10 1];ub=double(ub)';
lb = [0   0  0 0];lb=double(lb)';

myfun4 = @(b) (R-Rmax*(((x(:,1))+b(3))./(b(2)+(x(:,1))+b(1).*(x(:,2)))));
optnew=optimset(options,'MaxFunEvals',10000, 'MaxIter', 3000);
[b4,rss4] = lsqnonlin(myfun4,b0,lb,ub,optnew);
ICflag='AIC';

if j==1
    b1_best=b1;
    bic1 = biccustom(R,rss1,4,ICflag);
    min1=rss1;
    
    b2_best=b2;
    bic2 = biccustom(R,rss2,4,ICflag);
    min2=rss2;

    b3_best=b3;
    bic3 = biccustom(R,rss3,3,ICflag);
    min3=rss3;
    
    b4_best=b4;
    bic4 = biccustom(R,rss4,4,ICflag);
    min4=rss4;
end

if rss1<min1
    b1_best=b1;
    bic1 = biccustom(R,rss1,4,ICflag);
    min1=rss1;
end

if rss2<min2
    b2_best=b2;
    bic2 = biccustom(R,rss2,4,ICflag);
    min2=rss2;
end

if rss3<min3
    b3_best=b3;
    bic3 = biccustom(R,rss3,4,ICflag);
    min3=rss3;
end

if rss4<min4
    b4_best=b4;
    bic4 = biccustom(R,rss4,4,ICflag);
    min4=rss4;
end
end
AIC_sum=[bic1,bic2,bic3,bic4];
B.b1=b1;
B.b2=b2;
B.b3=b3;
b.b4=b4;