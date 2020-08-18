function AValues=loadAnalogChanel(filepath,StartTimes,chanNum)
Trialnumber = length(StartTimes)-1;
for i = 1:Trialnumber
currTrial = i;

    [adfreq, n, avalues] = plx_ad_span(filepath,chanNum-1,StartTimes(currTrial),StartTimes(currTrial+1)-1);
    currPhotoT = avalues;
    if length(currPhotoT) > 10000
          currPhotoT = currPhotoT(1:10000);
    end
    clear avalues adfreq n;
    if ~isempty(currPhotoT)
        AValues(currTrial,1:length(currPhotoT)) = currPhotoT';
    end
    
end
AValues(Trialnumber,10000)=0;