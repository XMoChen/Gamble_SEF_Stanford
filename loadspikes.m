function [ChannelNumber,NeuronNumber,tscounts]=loadspikes(filepath,StartTimes,outfileSpike)
[tscounts] = plx_info(filepath, 1);
tscounts(1,:) = [];
tscounts(:,1) = [];
ts1 = nonzeros(tscounts(:,1));
ts2 = nonzeros(tscounts(:,2));
ts3 = nonzeros(tscounts(:,3));
ts4 = nonzeros(tscounts(:,4));

% tscounts = nonzeros(tscounts);
% NeuronNumber = size(tscounts,1);
% ChannelNumber = size(tscounts,2);
l_ts1 = length(ts1);NeuronNumberCh1=l_ts1;
l_ts2 = length(ts2);NeuronNumberCh2=l_ts2;
l_ts3 = length(ts3);NeuronNumberCh3=l_ts3;
l_ts4 = length(ts4);NeuronNumberCh4=l_ts4;

    NeuronNumber = max([l_ts1,l_ts2,l_ts3,l_ts4]);

if ts4 ~= 0
    ChannelNumber = 4;
elseif ts3 ~= 0   
    ChannelNumber = 3;
elseif ts2 ~= 0   
    ChannelNumber = 2;
else
    ChannelNumber = 1;
end

Trialnumber = length(StartTimes)-1;
% load(outfile);
for currentchannel = 1:ChannelNumber
    for currentneuron = 0:NeuronNumber
        spikespertrial=[];currSpkT=[];
        %if ~isempty(tscounts(currentneuron,currentchannel))
        if (currentneuron==0 | tscounts(currentneuron,currentchannel) ~= 0)
            spike_ts=[];n=[];
            [n, spike_ts] = plx_ts(filepath, currentchannel, currentneuron);
            spike_ts = single(spike_ts*1000);
            for i = 1:Trialnumber
                currTrial = i;
                currSpkT = spike_ts(find((spike_ts >= StartTimes(currTrial)) ...
                        & (spike_ts <min([StartTimes(currTrial+1),StartTimes(currTrial)+3000]) )));
                    if ~isempty(currSpkT)
                        subtract = single(currSpkT - StartTimes(currTrial));
                        spikespertrial(currTrial,1:length(subtract)) = single(subtract)';
                        spikespertrial(find(spikespertrial < 0)) = 0;
                    end
            end
            
%             if exist(['SpikeTimes', num2str(currentneuron), num2str(currentchannel)],'var')
%             eval(['a=','SpikeTimes', num2str(currentneuron), num2str(currentchannel),';']);
%             [c1,l1]=size(a);[c2,l2]=size(spikespertrial);l=max(l1,l2);
%             if c1<Trialnumber
%             a(c1,l+1)=0;spikespertrial(c2,l+1)=0;
%             eval(['SpikeTimes', num2str(currentneuron), num2str(currentchannel), '= [a;spikespertrial];']);
%             else
%             eval(['SpikeTimes', num2str(currentneuron), num2str(currentchannel), '= a;'])
%             end
%             else
            eval(['SpikeTimes', num2str(currentneuron), num2str(currentchannel), '= single(spikespertrial);']); 
            save(outfileSpike,['SpikeTimes*'],'NeuronNumber','-append');
            clear SpikeTimes*
%             end
         end
    end
end

