function AlignEvent(file)
load(file);

EventName={'Target','Move','Result'};
times=round(Infortable{:,[24,25,27]});
I_Move=Infortable{:,19}~=0 & times(:,2)>400 & times(:,2)<10000 ;
I_Targ=Infortable{:,19}~=0 & times(:,1)>400 & times(:,2)<10000;
I_result=Infortable{:,19}~=0 & Infortable{:,20}~=0 & times(:,3)>400 & times(:,2)<10000;

N=size(Infortable,1);
Move=zeros(2,N,1024);
Target=zeros(2,N,1024);
Result=zeros(2,N,1024);
Fs=1000;
d = designfilt('bandstopiir','FilterOrder',2, ...
               'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
               'DesignMethod','butter','SampleRate',Fs);
LFPValues1 = filtfilt(d,LFPValues1')';
LFPValues2 = filtfilt(d,LFPValues2')';


for i=1:N
    
if I_Targ(i)==1  

Target(1,i,:)=LFPValues1(i,[(times(i,1)-500):(times(i,1)+523)]);
Target(2,i,:)=LFPValues2(i,[(times(i,1)-500):(times(i,1)+523)]);
end
if I_Move(i)==1
Move(1,i,:)=LFPValues1(i,[(times(i,2)-500):(times(i,2)+523)]);
Move(2,i,:)=LFPValues2(i,[(times(i,2)-500):(times(i,2)+523)]);
end
if I_result(i)==1
Result(1,i,:)=LFPValues1(i,[(times(i,3)-500):(times(i,3)+523)]);
Result(2,i,:)=LFPValues2(i,[(times(i,3)-500):(times(i,3)+523)]);
end
end

save(file,'Target','Move','Result','-append');
% plot(squeeze(mean(mean(Target,1),2)));
% plot(squeeze(mean(mean(Move,1),2)));
% plot(squeeze(mean(mean(Result,1),2)));