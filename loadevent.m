function [eventtimes,eventcodes,StartTimes]=loadevent(filepath)
[n, Time257, Code257] = plx_event_ts(filepath, 257);

event_ts = Time257 *1000; % to get ms from s
event_sv=Code257;

StartTimes = event_ts(find(Code257== 32512));
Trialnumber = length(StartTimes)-1;
loopNum=round(size(StartTimes ,1)/15);

    eventtimes= zeros(Trialnumber ,1);
    eventcodes= zeros(Trialnumber ,1);

for i = 1:Trialnumber

    currTrial = i;
    currEvT = event_ts(find((event_ts >= StartTimes(currTrial)) ...
            & (event_ts < StartTimes(currTrial+1))));
    currEvC = event_sv(find((event_ts >= StartTimes(currTrial)) ...
            & (event_ts < StartTimes(currTrial+1))));
   
            
       
        if ~isempty(currEvT)
            subtract = currEvT - StartTimes(currTrial);
            eventtimes(i,1:length(subtract)) = subtract';
            eventtimes(find(eventtimes< 0)) = 0;
        end
        if ~isempty(currEvC)
            eventcodes(i,1:length(currEvC)) = currEvC';
        end
        
        
        
end
eventtimes(size(eventtimes,1),24)=0;
eventcodes(size(eventtimes,1),24)=0;
