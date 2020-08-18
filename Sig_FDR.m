function Sig_I=Sig_FDR(Sig)
%size(Sig);
for i=1:size(Sig,1)
pID = FDR(Sig(i,:));
Sig_I(i,:)=Sig(i,:)<pID;
end