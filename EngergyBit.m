function B=EngergyBit(A)
% min(min(A)) 
% max(max(A))
% size(A)
% plot(abs(A));hold on;
% A=A+10^-3;
B=10*log10(abs(A));
B(B<-20)=-20;
