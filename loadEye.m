function [EyeX,EyeY]=loadEye(filepath,StartTimes)

EyeX=loadAnalogChanel(filepath,StartTimes,51);
EyeY=loadAnalogChanel(filepath,StartTimes,50);
% EyeV=sqrt(EyeX.^2+EyeY.^2);
% EyeV=medfilt1(EyeV,10);