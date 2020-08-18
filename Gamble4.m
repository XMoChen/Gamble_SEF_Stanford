clear all; close all;
%%%%%%%%%%%%%%%%%%%%
%table  1. display pattern   2. max reward1     3. max reward2
%       4. trial type        5. prob1           6. prob 2   
%       7. Atractor 1        8. Atractor 2      9. Image pattern 1
%       10 Image pattern 2   11. Fixation page On  12. Fix  In
%       13 Atractor page On  14. Target Page On    15. Target In
%       16.Choice            17. Chosen Target On  18. Gamble Result
%       19. Result On        20. Reward On         21. Reward Amount
%       22. Image1           23. Image2            24. Chosen Image
%       25. NonChosen Image  26. Move End          27 Move direction
%       28. Chosen Onset     29. nonChosen Onset   30. SubV 1
%       31. SubjV 2          32. Chosen SubjV      33.nonChosen SubjV
%       34. Rt               35, temperature for cooling
%%%%%%%%%%%%%%%%%%%%%%
 
%1 gamble task without cooling
%  FlagName={'go09022011(2)','go09022011(3)','go09052011','go09072011(2)','go09112011','go09122011','go09142011','go09162011',...  
%     'go09192011','go09212011-01','go09232011','go09272011-01','go09292011-01','go09302011-01','go10022011_mrg-02','go10032011-01','go10042011-01','go10052011','go10062011-02','go10102011-01','go10112011-01','go10122011-01','go01172012-02',...
%     'go10192011-01','go10202011-01','go01192012-02',...
%     'go01232012-01','go01262012-01','go01272012-01','go01292012-01','go01302012-01','go02052012-01'};
%  FlagName={'cgonew04242012-01_sub1','cgonew04242012-01_sub2','cgonew05102012-01_sub1','cgonew05102012-01_sub2','Ago03082013Newb'};

% % single channel 
% multichannels 
% 8 positions   train 513

%2** cooling experiment try
%  FlagName={'cgo02242012_mrg-01','cgo02272012_mrg','cgo02272012c_mrg','cgo03012012_mrg','cgo03042012_mrg','cgo03062012'};
%2 cooling experiment 8 targets
%      FlagName={'cgo03082012-01','cgo03112012-01','cgo03132012-01','cgo03152012-01','cgo03202012-02','cgo03252012-01','cgo03272012-01','cgo03292012-01','cgo04012012-01','cgo04032012-01'};
  %%%pool short long together
%   FlagName={'cgo03082012-01','cgo03112012-01','cgo03132012-01','cgo03152012-01','cgo03202012-02','cgo03252012-01','cgo03272012-01','cgo03292012-01','cgo04012012-01','cgo04032012-01',...
%             'cgolong04122012t-01','cgolong04152012-01','cgolong04172012-01','cgo04222012long-01',...
 FlagName={ 'cgonew04242012-01','cgonew05102012-01','cgonew05122012','cgonew05132012-01','cgonew05152012-01'};

% FlagName={ 'I10222012cool','Icool10262012-01','Igocooling11022012_mrg','Igocooling11022012c','Igocool11122012-01','Igocool11142012-01','Igocool11142012b-01','Igocool11182012-01','Igocool11182012b','Igocool12032012','Icgo12052012new','Icgo12052012new(b)','Icgo12102012new_mrg','Icgo12122012new','Icgo12122012new(b)'};

        
        %   depth=[111,718;442,1343;94,832;142,-113;618,525;778,-345;600,1500;1075,575;772,983;-44,686;...
%          0,320;392,-33;732,637;823,1023;...
%          410,71;483,1132;250,75;688,254;253,382];
%   depthout=[-740,34;0,0;-230,-57;-1030,-900;-466,-600;-100,-940;0,0;137,-940;-200,-300;-1100,-300;...
%          -634,-630;-787,-670;-200,-700; -330,-700;...
%          -300,-200;0,400;-550,-730;-199,-200;-80,-300];
%  depth=depth-depthout; [depthS,depthI]=sort(depth);
%   A=reshape(depth',1,38);
     
%3 test the subjective value of sure red-3 and sure white-4
%   FlagName={'gored04102012','gowhite04102012'};  
% test cooling %Flag={'cgolong04102012'};

% 4 document the effect of cooling
%     FlagName={'cgolong04122012t-01','cgolong04152012-01','cgolong04172012-01','cgo04222012long-01'};

% 5 document the subjective value of the new targets  color reversed
%   FlagName={'cgonew04242012-01','cgonew05102012-01','cgonew05122012','cgonew05132012-01','cgonew05152012-01'};

% 6 learn new sure targets  
% FlagName={'Sure605172012','Sure605172012(2)','Sure605172012(3)','Sure605182012(1)','Sure605182012(2)','Sure605182012(3)','Sure605212012','Sure605232012'};
% FlagName={ 'SG605212012','SG605232012','GS05242012'};


mainpath='D:\Projects\SEFCooling\data\';
Plexonpath='D:\Projects\SEFCooling\Aragorn plexon data\';
% outfile=['I:\data\',char(Flag(1))];



eventtimes=[];eventcodes=[];LFPValues=[];
EyeV=[];NeronAll=0;StempAll=[];

for i=1:size(FlagName,2)%23
Flag0=char(FlagName(i))
% Creatfolders(Flag0,mainpath);
outfile=[mainpath,'data\',char(FlagName(i)),'\',char(FlagName(i)),'.mat'];
outfileSpike=[mainpath,'data\',char(FlagName(i)),'\',char(FlagName(i)),'Spike.mat'];
outfileLFP=[mainpath,'data\',char(FlagName(i)),'\',char(FlagName(i)),'LFP.mat'];
figurepath=[mainpath,'figure\',Flag0,'\'];
filepath=[Plexonpath,char(FlagName(i)),'.plx'];  
% 
% LFPIn=[mainpath,'LFP1\in\',Flag0,'\'];
% LFPOut=[mainpath,'LFP1\out\',Flag0,'\'];
% LFPtrial=[mainpath,'LFP1\trial\',Flag0,'\'];
% if exist(LFPIn)==0
% mkdir(LFPIn);mkdir(LFPOut);mkdir(LFPtrial);
% end
% 
% % %%%%%%%%%%%%%%%%%%%
%  [n, Time257, Code257] = plx_event_ts(filepath, 257);
% 
% [eventtimes,eventcodes,StartTimes]=loadevent(filepath);
% save(outfile,'eventtimes','eventcodes','StartTimes','-append');
% % % % 
%   load(outfile);
% % %%%%Spike 
%  [ChannelNumber,NeuronNumber,tscounts]=loadspikes(filepath,StartTimes,outfileSpike);
% %%%LFP 
%  [ChannelNumber,NeuronNumber]=channelinfo(filepath,outfile);
% for ii=1:ChannelNumber
% eval(['LFPValues',num2str(ii),'=loadAnalogChanel(filepath,StartTimes,',num2str(ii),');']);
% save(outfileLFP,['LFPValues*'],'-append');
% clear LFPValues*;
% end
% % %%% saccade signal
load(outfile,'eventtimes','eventcodes','StartTimes');
[EyeX,EyeY]=loadEye(filepath,StartTimes);
save(outfileLFP,'EyeX','EyeY','-append');
clear Eye*;
  load(outfile);
% % %%%%%%%%%%%%%%%%%%   
% 
% 

% % TrialNum=size(eventtimes,1);
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %1  for recording experiemnt
% imageid1=[4,6,7,3, 4,6,7,3, 1:7       2:7         3:7         4:7          5:7         6:7            7    8 ];
% imageid2=[2,5,3,1, 2,5,3,1, zeros(1,7)   ones(1,6)   2*ones(1,5) 3*ones(1,4)  4*ones(1,3) 5*ones(1,2) 6    0 ];

% table=EventOrganize(imageid1,imageid2,eventcodes,eventtimes,outfile);
%  MoveTime(outfileLFP,outfile);
%  BehaviorAnalysis(imageid1,imageid2,table,figurepath,outfile);
%  SpikeDensity(outfile,'Targ',400,100,outfileSpike);
%  SpikeDensity(outfile,'Move',400,200,outfileSpike);
% SpikeDensity(outfile,'BackGround',400,400,outfileSpike);
%  SpikeDensity(outfile,'Reward',500,200,outfileSpike);
% 
% 
% %  NeronN=SpikeDirection2(outfile,'Move',400,200,outfileSpike,figurepath,'NoChoice');
% % NeronAll=NeronAll+NeronN;
% SpikeDirectionValue0(outfile,'Move',400,200,outfileSpike,figurepath,'NoChoice');
%  SpikeDirectionValue0(outfile,'Move',400,200,outfileSpike,figurepath,'Choice');
% SpikeDirectionDV(outfile,'Move',400,200,outfileSpike,figurepath,'Choice');


% SpikeDirectionValueR(outfile,'Reward',500,200,outfileSpike,figurepath,'NoChoice');
% SpikeDirectionValueR(outfile,'Reward',500,200,outfileSpike,figurepath,'Choice');

%  Alignflag='Move';length=1400;off=700;
%  LFPdataOrganize(outfile,outfileLFP,LFPtrial,Alignflag,length,off,LFPIn)
%  Alignflag='Targ';length=1400;off=700;
%  LFPdataOrganize(outfile,outfileLFP,LFPtrial,Alignflag,length,off,LFPIn)
%

%2 for early colling experiment 4 targets without recording temperature
% 'cgo02242012_mrg-01','cgo02272012_mrg','cgo02272012c_mrg'
% imageid1=[ 5:7  6:7  7  5:7  6:7  7  5:7  6:7  7  5:7  6:7  7  4 5 6 7    9 ];
% imageid2=[ 4*ones(1,3) 5*ones(1,2) 6 4*ones(1,3) 5*ones(1,2) 6 4*ones(1,3) 5*ones(1,2) 6 4*ones(1,3) 5*ones(1,2) 6 0 0 0 0    0 ];
%  table=EventOrganizeCooling2(imageid1,imageid2,eventcodes,eventtimes,outf
%  ile);  % without temperatuer

%3 for early colling experiment  4 targets
%'cgo03012012_mrg' 'cgo03042012_mrg' 'cgo03062012_mrg'
% imageid1=[ 5:7  6:7  7  5:7  6:7  7  5:7  6:7  7  5:7  6:7  7  4 5 6 7   9  8  ];
% imageid2=[ 4*ones(1,3) 5*ones(1,2) 6 4*ones(1,3) 5*ones(1,2) 6 4*ones(1,3) 5*ones(1,2) 6 4*ones(1,3) 5*ones(1,2) 6 0 0 0 0  0  0 ];
%   table=EventOrganizeCooling4(imageid1,imageid2,eventcodes,eventtimes,out
%   file);  % with temperature

%4 for late cooing experiment  7 targets   with sure target block
% imageid1=[4,6,7,3, 4,6,7,3, 1:7       2:7         3:7         4:7          5:7         6:7            7    8  9];
% imageid2=[2,5,3,1, 2,5,3,1, zeros(1,7)   ones(1,6)   2*ones(1,5) 3*ones(1,4)  4*ones(1,3) 5*ones(1,2) 6    0  0];
% BehaviorAnalysisCoolingProb(imageid1,imageid2,figurepath,outfile,Flag0)
% table=EventOrganizeCooling4(imageid1,imageid2,eventcodes,eventtimes,outfile);
% MoveTime(outfileLFP,outfile);
% BehaviorAnalysisCooling(imageid1,imageid2,table,outfile); 
% ChoiceProbCool2(outfile);
%    SpikeDensity(outfile,'Targ',500,200,outfileSpike);
% SpikeDensity(outfile,'Move',400,200,outfileSpike);
% SpikeDensity(outfile,'Reward',500,200,outfileSpike);
%  [PsiV_Comp,PsiV_SD_Comp]=BehaviorCoolingSubV(imageid1,imageid2,outfile)
%   [Stemp]=SpikeAllCooling(outfile,outfileSpike,depth(i,:),figurepath);

%  NeronN=SpikeDirectionCool2(outfile,'Move',400,200,outfileSpike,figurepath,'NoChoice');
% SpikeDirectionValueCool0(outfile,'Move',400,200,outfileSpike,figurepath,'NoChoice');
%  SpikeDirectionValueCool0(outfile,'Move',400,200,outfileSpike,figurepath,'Choice');

% SpikeDirectionValueCoolR(outfile,'Reward',400,200,outfileSpike,figurepath,'NoChoice');
% SpikeDirectionValueCoolR(outfile,'Reward',400,200,outfileSpike,figurepath,'Choice');

%  Alignflag='Move';length=1400;off=700;
%  LFPdataOrganize(outfile,outfileLFP,LFPtrial,Alignflag,length,off,LFPIn)
%  Alignflag='Targ';length=1400;off=700;
%  LFPdataOrganize(outfile,outfileLFP,LFPtrial,Alignflag,length,off,LFPIn)
%

%    LFPBaseCooling(outfile,LFPtrial,outfileLFP)
%   LFPAllCooling(outfile,LFPtrial,outfileSpike,depth,figurepath,outfileLFP)
%  figurepath='D:\Gamble\data\figure\coolingLFP\';   
% LFPAllCoolingT(outfile,LFPtrial,outfileSpike,depth,figurepath,outfileLFP,Flag0)
%   LFPAllCoolingTTemp(outfile,LFPtrial,outfileSpike,depth,figurepath,outfileLFP,Flag0)
 
% SpikeMulti_D_V(outfileSpike,outfile);
% LFPDirection(outfile,LFPtrial,outfileSpike,depth,figurepath,outfileLFP,Flag0)

  
  
% SpikeDensity(outfile,'BackGround',400,400,outfileSpike);
%      Stemp=SpikeAllCooling(outfile,outfileSpike,depth(i,:));
%     StempAll=[StempAll;Stemp];
%     StempAllS=StempAll(depthI,:);
%     for j=1:3
%        clear I
%        I=(A<600*j & A>600*(j-1)) ;
%        StempAllS2(j,:)=mean(StempAll(I',:));
%     end

% to document the effect of cooling in long interval
% BehaviorAnalysisCoolingE22long.m

 %5 for subjective value of white (1) and red (0.44) sure target
% imageid1=[1:7  8       2:7         3:7         4:7          5:7         6:7            7    1:7       9  8];
% imageid2=[zeros(1,7) 0   ones(1,6)   2*ones(1,5) 3*ones(1,4)  4*ones(1,3) 5*ones(1,2) 6   8*ones(1,7) 0  0];
% table=EventOrganizeCooling4(imageid1,imageid2,eventcodes,eventtimes,outfile);
% [PsiV,ProT]=BehaviorAnalysis8(imageid1,imageid2,table,figurepath)

%6 for cooing experiment 8   choice new targets
imageid1=[4,6,7,3, 4,6,7,3, 1:7       2:7         3:7         4:7          5:7         6:7            7    9  8];
imageid2=[2,5,3,1, 2,5,3,1, zeros(1,7)   ones(1,6)   2*ones(1,5) 3*ones(1,4)  4*ones(1,3) 5*ones(1,2) 6    0  0];

%7 for comparison of the sure targets 1,2,3,4,5,9
% imageid1= [1,2,3,4,5,6,   2,3,4,5,6,   3,4,5,6,  4,5,6,   5,6,   6];
% imageid2= [0,0,0,0,0,0,         1,1,1,1,1,   2,2,2,2,  3,3,3,   4,4,   5];


% for Sure Gamble comparison
% imageid1= [8,9,10,11,12,13,  1,2,3,4,5,6,7,   8,9,10,  8,9,10,  8,9,10,11,12,   8,9,10,11,12,   8,9,10,11,12,  8,10,12,13,  8,10,12,13];
% imageid2= [0,0,0,0,0,0,        0,0,0,0,0,0,0,   1,1,1,     2,2,2,      3,3,3,3,3,        4,4,4,4,4,        5,5,5,5,5,       6,6,6,6,      7,7,7,7];


% table=EventOrganizeCooling4(imageid1,imageid2,eventcodes,eventtimes,outfile);
% load(outfile)
% BehaviorAnalysis(imageid1,imageid2,table,figurepath);
% PsiVC=BehaviorAnalysisCoolingE2(imageid1,imageid2,table,figurepath,'Cool');
% PsiVN=BehaviorAnalysisCoolingE2(imageid1,imageid2,table,figurepath,'Normal');

% % Movement time: MoveOn MoveEnd
% table=MoveTime(outfileLFP,table,outfile);


% PsiVC=BehaviorAnalysisCoolingE2(imageid1,imageid2,table,figurepath,'Cool');
% PsiVN=BehaviorAnalysisCoolingE2(imageid1,imageid2,table,figurepath,'Normal');
% save(outfile,'PsiVC','PsiVN','-append');

% PsiV=BehaviorAnalysis(imageid1,imageid2,table,figurepath);
% PsiV=BehaviorAnalysisCooling(imageid1,imageid2,table,figurepath);
% 
% save(outfile,'imageid1','imageid1','eventtimes','eventcodes','TrialNum','PsiV','table','-append');
% infotable=FunInfotable(outfile);
% clear 'eventtimes'
%  SpikeDirectionCooling(outfile,'Targ',outfileSpike,figurepath,'TwoChoice','Cool');  % Move Targ
% SpikeDirection(outfile,'Targ',outfileSpike,figurepath,'NoChoice');  % Move Ta

 % the other direction is fixed, therefore correlated with chosen direction
% SpikeNChDirection(outfile,'Choice',figurepath);  

% trialflag='NoChoice'; % only no choice trial  1 no choice 2 two choice
% SpikeValue(outfile,'Move',figurepath,outfileSpike,trialflag);   % Move Targ Choice
% % SpikeNChValue(outfile,'Choice',figurepath,outfileSpike,trialflag);   % Value of the nonchosen target Move Targ Choice
% SpikeValue(outfile,'Result',figurepath,outfileSpike,trialflag);   % Move Targ Choice
% 
% 
% SpikeResult(outfile,'Result',figurepath,outfileSpike)   % Result  Reward

% SpikeValueDirection(outfile,'Choice',figurepath)   % Targ Move Choice


% begin to look at the onset effect
% flagtrial=1;
% SpikeOnset(outfile,'Targ',figurepath,outfileSpike,flagtrial)
% SpikeOnset(outfile,'Move',figurepath,outfileSpike,flagtrial)
 end
% 
% 
% % sure vs gamble 1 lose and l sure
