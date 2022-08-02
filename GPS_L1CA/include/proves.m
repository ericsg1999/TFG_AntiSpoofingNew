%[data,fid] = readSignalFile(fileID,pos,settings,numSamples,openFile)settings = initSettings();
APT_period=6; %mostres

settings=initSettings();
[data1,fid] = readSignalFile(0,0,settings,3,1)
ftell1=ftell(fid)

%pos=0.00000004*settings.samplingFreq;
it=1;
disp('----------------------------------')
[data2,fid] = readSignalFile(fid,(it*APT_period-1)*4,settings,3,0)
ftell2=ftell(fid)

it=2;
disp('----------------------------------')
[data3,fid] = readSignalFile(fid,(it*APT_period-1)*4,settings,3,0)
ftell3=ftell(fid)
%% 
settings=initSettings();
samples=1050100;
[data1,fid] = readSignalFile(0,0,settings,samples,1);
a=2;
%%
a=[1 2 3; 4 5 6; 7 8 9]
a(:,(1:2))