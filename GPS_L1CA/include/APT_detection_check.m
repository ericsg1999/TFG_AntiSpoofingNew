function [fileID,spoofingRecord,secondaryPeakMagnitudeVector] = APT_detection_check(settings,fid1,filePos,SatellitePresentList,openFile,spoofingRecord,secondaryPeakMagnitudeVector)
%APT_DETECTION_CHECK calls acquisition function in order to perform the APT
%check with the purpose of detecting a possible secondary peak in the
%search grid. 


%% Read all the samples in an APT period
%-----------Compute samples in an APT period-------------------------------
samplesPerCode = round(settings.samplingFreq / ...
                           (settings.codeFreqBasis / settings.codeLength));

%numSamples=settings.AptPeriod*samplesPerCode;
codeNum=42;
numSamples=codeNum*samplesPerCode; 
% % % % %----------Read Signal File------------------------------------------------
% % % % %raw_signal_AptPeriod_long = fread(fid1, settings.AptPeriod*samplesPerCode, settings.dataType)'; %settings.AptPeriod corresponds to the ms between apt check and samplesPercode, as 1 code=1ms, correspond to the samples per code. the multiplication is the samples between apt check
% % % % [raw_signal_AptPeriod_long fileID]=readSignalFile(fid1,settings,numSamples,openFile);
% % % % %ftell(fileID)
% % % % raw_signal_42ms=raw_signal_AptPeriod_long(1:42*samplesPerCode);

%----------Read Signal File------------------------------------------------
%raw_signal_AptPeriod_long = fread(fid1, settings.AptPeriod*samplesPerCode, settings.dataType)'; %settings.AptPeriod corresponds to the ms between apt check and samplesPercode, as 1 code=1ms, correspond to the samples per code. the multiplication is the samples between apt check
[raw_signal_AptPeriod_long fileID]=readSignalFile(fid1,filePos,settings,numSamples,openFile);
%ftell(fileID)
raw_signal_42ms=raw_signal_AptPeriod_long(1:42*samplesPerCode);

%% Perform APT acquisition
acqType='APT';
[acqResults,spoofingAlert,secondaryPeakMagnitude] = acquisition_new_w_APT(raw_signal_42ms, settings, acqType,SatellitePresentList);

spoofingRecord=[spoofingRecord spoofingAlert];
secondaryPeakMagnitudeVector=[secondaryPeakMagnitudeVector secondaryPeakMagnitude];
end

%[channel.PRN]