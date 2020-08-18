%
close all;
x=0:0.01:2*pi;
v=ones(size(x));
y=0.035+ sigmoid2([28,0.59,0.1],v).*gaussian2([1,0.8,1.6],x);
plot(x,y); hold on
y=0.035+ sigmoid2([28,0.59,0.1],0.5*v).*gaussian2([1,0.8,1.6],x);
plot(x,y,'g'); hold on
y=0.035+ sigmoid2([28,0.59,0.1],0.3*v).*gaussian2([1,0.8,1.6],x);
plot(x,y,'k'); hold on