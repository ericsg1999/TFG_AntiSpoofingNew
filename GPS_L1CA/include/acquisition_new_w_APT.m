function acqResults = acquisition_new_w_APT(longSignal, settings, acqType,SatellitePresentList)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%       longSignal    - 11 ms of raw signal from the front-end
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number.

%--------------------------------------------------------------------------
%                         CU Multi-GNSS SDR
% (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
% Based on the original work by Darius Plausinaitis,Peter Rinder,
% Nicolaj Bertelsen and Dennis M. Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

%CVS record:
%$Id: acquisition.m,v 1.1.2.12 2006/08/14 12:08:03 dpl Exp $

%% Condition input signal to speed up acquisition ===================

% If input IF signal freq. is too high, a resampling strategy is applied
% to speed up the acquisition, whcih is selectable.
if (settings.samplingFreq > settings.resamplingThreshold && ...
        settings.resamplingflag == 1)
    
    %--- Filiter out signal power outside the main lobe of CM code --------
    fs = settings.samplingFreq;
    IF = settings.IF;
    % Bandwidth of CA mian lobe: 0.5e6 is a margin to make sure most of CA
    % code power will not filtered out
    BW = settings.codeFreqBasis*2 + 0.5e6;
    % Filter parameter
    w1 = (IF)-BW/2;
    w2 = (IF)+BW/2;
    wp = [w1*2/fs w2*2/fs];
    % Filter coefficients
    b  = fir1(700,wp);
    % Filter operation
    longSignal = filtfilt(b,1,longSignal);
    
    % --- Find resample frequency -----------------------------------------
    % Refer to bandpass sampling theorem(Yi-Ran Sun,Generalized Bandpass
    % Sampling Receivers for Software Defined Radio)
    
    % Upper boundary frequency of the bandpass IF signal
    fu = settings.IF + BW/2;
    % Lower freq. of the acceptable sampling Freq. range
    n = floor(fu/BW);
    if (n<1)
        n = 1;
    end
    lowerFreq = 2*fu/n;
    % Lower boundary frequency of the bandpass IF signal
    fl = settings.IF - BW/2;
    % Upper boundary frequency of the acceptable sampling Freq. range
    if(n>1)
        upperFreq = 2*fl/(n-1);
    else
        upperFreq = lowerFreq;
    end
    % Save orignal Freq. for later use
    oldFreq = settings.samplingFreq;
    
    % Take the center of the acceptable sampling Freq. range as
    % resampling frequency. As settings are used to generate local
    % CM code samples, so assign the resampling freq. to settings.
    % This can not change the settings.samplingFreq outside this
    % acquisition function.
    settings.samplingFreq = ceil((lowerFreq + upperFreq)/2);
    
    %--- Downsample input IF signal ---------------------------------------
    % Signal length after resampling
    signalLen = floor((length(longSignal)-1) /oldFreq * settings.samplingFreq);
    % Resampled signal index
    index = ceil((0:signalLen-1)/settings.samplingFreq * oldFreq);
    index(1) = 1;
    % Resampled signal
    longSignal = longSignal(index);
    % For later use
    oldIF = settings.IF;
    % Resampling is equivalent to down-converting the original IF by integer
    % times of resampling freq.. So the IF after resampling is equivalent to:
    settings.IF = rem(settings.IF,settings.samplingFreq);
end % resampling input IF signals

%% Initialization ===================================================
%--- Varaibles for coarse acquisition -------------------------------------
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
% Find sampling period
ts = 1 / settings.samplingFreq;
% Find phase points of 2ms local carrier wave (1ms for local duplicate,
% the other 1ms for zero padding)
phasePoints = (0 : (samplesPerCode * 2 -1)) * 2 * pi * ts;
% Number of the frequency bins for the specified search band
numberOfFreqBins = round(settings.acqSearchBand * 2 / settings.acqSearchStep) + 1;
% Carrier frequency bins to be searched
coarseFreqBin = zeros(1, numberOfFreqBins);

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
switch acqType
    case 'Normal'
        acqResults.carrFreq     = zeros(1, 32);
        % C/A code phases of detected signals
        acqResults.codePhase    = zeros(1, 32);
        % Correlation peak ratios of the detected signals
        acqResults.peakMetric   = zeros(1, 32);
    case 'APT'
        acqResults.carrFreq     = zeros(settings.AptNumberChannelsPerSat, 32);
        % C/A code phases of detected signals
        acqResults.codePhase    = zeros(settings.AptNumberChannelsPerSat, 32);
        % Correlation peak ratios of the detected signals
        acqResults.peakMetric   = zeros(settings.AptNumberChannelsPerSat, 32);
end

acqResults.SatellitePresentList= [];

%--- Varaibles for fine acquisition ---------------------------------------
% Carrier frequency search step for fine acquisition
fineSearchStep = 25;
% Number of the frequency bins for fine acquisition
numOfFineBins = round(settings.acqSearchStep/ fineSearchStep) + 1;
% Carrier frequencies of the fine frequency bins
fineFreqBins = zeros(1, numOfFineBins);
% Correlation values for all fine frequency bins
fineResult = zeros(1,numOfFineBins);
% Coherent integration for each of 40 codes
sumPerCode = zeros(1,40);
% Phase points of the local carrier wave
finePhasePoints = (0 : (40*samplesPerCode-1)) * 2 * pi * ts;

%--- Input signal power for GLRT statistic calculation --------------------
sigPower = sqrt(var(longSignal(1:samplesPerCode)) * samplesPerCode);

% Perform search for all listed PRN numbers ...
%fprintf('(');
list_sat_to_acquire=zeros(1,32);
switch acqType
    case 'Normal'
        list_sat_to_acquire=settings.acqSatelliteList;
    case 'APT'
        list_sat_to_acquire=SatellitePresentList;
    otherwise 
        fprintf('ERROR: Incorrect Acquisition type (acqType) input parameter of the function acquisition.m');
end
    
for PRN = list_sat_to_acquire
    %% Coarse acquisition ===========================================
    % Generate C/A codes and sample them according to the sampling freq.
    caCodesTable = makeCaTable(PRN,settings);
    % Add zero-padding samples
    caCodes2ms = [caCodesTable zeros(1,samplesPerCode)];
    % Search results of all frequency bins and code shifts (for one satellite)
    results = zeros(numberOfFreqBins, samplesPerCode*2);
    %--- Perform DFT of C/A code ------------------------------------------
    caCodeFreqDom = conj(fft(caCodes2ms));
    
    %--- Make the correlation for all frequency bins
    for freqBinIndex = 1:numberOfFreqBins
        % Generate carrier wave frequency grid
        coarseFreqBin(freqBinIndex) = settings.IF + settings.acqSearchBand - ...
            settings.acqSearchStep * (freqBinIndex - 1);
        % Generate local sine and cosine
        sigCarr = exp(-1i * coarseFreqBin(freqBinIndex) * phasePoints);
        
        %--- Do correlation -----------------------------------------------
        for nonCohIndex = 1: settings.acqNonCohTime
            % Take 2ms vectors of input data to do correlation
            signal = longSignal((nonCohIndex - 1) * samplesPerCode + ...
                1 : (nonCohIndex + 1) * samplesPerCode);
            % "Remove carrier" from the signal
            I      = real(sigCarr .* signal);
            Q      = imag(sigCarr .* signal);
            % Convert the baseband signal to frequency domain
            IQfreqDom = fft(I + 1i*Q);
            % Multiplication in the frequency domain (correlation in
            % time domain)
            convCodeIQ = IQfreqDom .* caCodeFreqDom;
            % Perform inverse DFT and store correlation results
            cohRresult = abs(ifft(convCodeIQ));
            % Non-coherent integration
            results(freqBinIndex, :) = results(freqBinIndex, :) + cohRresult;
        end % nonCohIndex = 1: settings.acqNonCohTime
    end % frqBinIndex = 1:numberOfFreqBins
    
    %% Look for correlation peaks for coarse acquisition ============
%     % Find the correlation peak and the carrier frequency
%     [~, acqCoarseBin] = max(max(results, [], 2));
%     % Find code phase of the same correlation peak
%     [peakSize, codePhase] = max(max(results));
%      
%     
%     
    %--------------- FIRST PEAK -------------------------------------------
    %--- Find the correlation peak and the carrier frequency --------------
    %M = max(A,[],dim) returns the maximum element along dimension dim. For example, if A is a matrix, then max(A,[],2) is a column vector containing the maximum value of each row.
    [primaryPeakSize frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [primaryPeakSize primaryCodePhase] = max(max(results));
    
    codephase_inferiorLimit=1;
    codephase_superiorLimit=size(results,2);
    
    %if primaryCodePhase is in the 1st ms, look for the second and third peak in
    %the 1st ms. If primaryCodePhase is in the 2nd ms, look for the
    %second and third peak in the 2nd ms
    if primaryCodePhase<=samplesPerCode %the first peak is in the 1st ms
        codephase_inferiorLimit=1;
        codephase_superiorLimit=samplesPerCode;
    else %the first peak is in the 2nd ms
        codephase_inferiorLimit=samplesPerCode;
        codephase_superiorLimit=size(results,2);
    end
    
    %--------------- SECOND PEAK ------------------------------------------
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
%     samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    samplesPerCodeChip=round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = primaryCodePhase - samplesPerCodeChip;
    excludeRangeIndex2 = primaryCodePhase + samplesPerCodeChip;
    
    
    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < codephase_inferiorLimit+2
        codePhaseRange = excludeRangeIndex2 : ...
                         (codephase_superiorLimit-1 + excludeRangeIndex1);
                         
    elseif excludeRangeIndex2 >= codephase_superiorLimit
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [codephase_inferiorLimit:excludeRangeIndex1, ...
                          excludeRangeIndex2 : codephase_superiorLimit];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    [secondaryPeakSize junk] = max(results(frequencyBinIndex, codePhaseRange));
    
    %Find the codephase of the second peak outside the surroundings of the
    secondaryPeakCodePhase=find(results(frequencyBinIndex,:)==secondaryPeakSize);
    
    %--------------- THIRD PEAK -------------------------------------------
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    excludeRangeIndex11 = secondaryPeakCodePhase - samplesPerCodeChip;
    excludeRangeIndex22 = secondaryPeakCodePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange2 = [excludeRangeIndex2 : excludeRangeIndex11...
                          excludeRangeIndex22:(codephase_superiorLimit-1 + excludeRangeIndex1)];

    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange2 = [(excludeRangeIndex2 - samplesPerCode):excludeRangeIndex11 ...
                         excludeRangeIndex22:excludeRangeIndex1] ;
    else
        if excludeRangeIndex1<excludeRangeIndex11
            codePhaseRange2 = [codephase_inferiorLimit:excludeRangeIndex1  excludeRangeIndex2:excludeRangeIndex11 excludeRangeIndex22:codephase_superiorLimit];
        else
            codePhaseRange2 = [codephase_inferiorLimit:excludeRangeIndex11  excludeRangeIndex22:excludeRangeIndex1 excludeRangeIndex2:codephase_superiorLimit];
        end
    end
    
    %--- Find the third highest correlation peak in the same freq. bin ---
    [thirdPeakSize junk] = max(results(frequencyBinIndex, codePhaseRange2));
    
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(1,PRN) = primaryPeakSize/thirdPeakSize;
    acqResults.codePhase(1,PRN) = primaryCodePhase;
    %acqResults.peakMetric(1,PRN) = peakSize/sigPower/settings.acqNonCohTime;
    
    
    
    %% plot the grid of the acquired acquisition search grids
 
    if settings.acqInitialPlots==1 && (strcmp(acqType,'Normal')==1)
        Td=0:1:(samplesPerCode-1);
        figure;
        if primaryCodePhase<samplesPerCode
            mesh(Td,coarseFreqBin,results(:,(1:samplesPerCode)));
        else
            mesh(Td,coarseFreqBin,results(:,(samplesPerCode+1:2*samplesPerCode)));
        end
        title(['PCPS ' acqType ' Acquisition grid for SV ID ', num2str(PRN)]);
        xlabel('Code delay [samples]');
        ylabel('Doppler freq [Hz]');
        
    elseif settings.AptPlots==1 && (strcmp(acqType,'APT')==1)
        Td=0:1:(samplesPerCode-1);
        
        if settings.AptShowPlots==0
            f=figure;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
            
%             f=figure
%             frame_h = get(handle(gcf),'JavaFrame');
%             set(frame_h,'Maximized',1);
            set(f, 'visible', 'off');
            
        else
            f=figure
        end
        if primaryCodePhase<samplesPerCode
            mesh(Td,coarseFreqBin,results(:,(1:samplesPerCode)));
        else
            mesh(Td,coarseFreqBin,results(:,(samplesPerCode+1:2*samplesPerCode)));
        end
        titol=['PCPS ' acqType ' Acquisition grid for SV ID ', num2str(PRN) ' (', datestr(datetime('now')) ')'];
%         data_string=strsplit(datestr(datetime('now')))
%         titol_=['PCPS_' acqType '_PRN_', num2str(PRN) '_(', cell2mat(split(datestr(datetime('now')),' ',2)) ')'];
        titol_=['PCPS_' acqType '_PRN_', num2str(PRN) '_(', strrep(strrep(datestr(now), ':', '-'),' ','-'), ')'];
        title(titol);
        xlabel('Code delay [samples]');
        ylabel('Doppler freq [Hz]');
        
        if settings.AptSavePlots==1
            %set(f,'WindowState','fullscreen');
            saveas(f,['D:\Images\' titol_ '.fig']); %gcf stands for the actual figure
        end
        
    end
    
    
    
    
    
    
    % If the result is above threshold, then there is a signal ...
    %% Fine carrier frequency search ================================
    
    %--- Do fine acquisition -----------------------------------
    if acqResults.peakMetric(1,PRN) > settings.acqThreshold
        
        % Indicate PRN number of the detected signal
        %fprintf('%02d', PRN);
        
        %--- Prepare 20ms code, carrier and input signals -----------------
        % C/A code with 10230 chips
        caCode = generateCAcode(PRN);
        % C/A code sample index
        codeValueIndex = floor((ts * (0 : 40*samplesPerCode -1)) / ...
            (1/settings.codeFreqBasis));
        % C/A code samples
        caCode40ms = caCode(rem(codeValueIndex, settings.codeLength) + 1);
        
        % Take 40cm incoming signal for fine acquisition
        sig40cm = longSignal(primaryCodePhase:primaryCodePhase + 40*samplesPerCode -1);
        
        %--- Search different fine freq bins ------------------------------
        for fineBinIndex = 1 : numOfFineBins
            %--- Correlation for each code --------------------------------
            % Carrier frequencies of the frequency bins
            fineFreqBins(fineBinIndex) = coarseFreqBin(frequencyBinIndex) + ...
                settings.acqSearchStep/2 - fineSearchStep * (fineBinIndex - 1);
            % Local carrier signal
            sigCarr40cm = exp(-1i * fineFreqBins(fineBinIndex) * finePhasePoints);
            % Wipe off code and carrier from incoming signals
            basebandSig = sig40cm .* caCode40ms .* sigCarr40cm;
            
            % Coherent integration for each code
            for index = 1:40
                sumPerCode(index) = sum(basebandSig( samplesPerCode * ...
                    (index - 1) + 1 : samplesPerCode * index ));
            end
            
            %--- Search Nav bit edge for ----------------------------------
            % 20 cases of Nav bit edge
            maxPower = 0;
            for comIndex = 1:20
                % Power for 20ms coherent integration
                comPower = abs(sum(sumPerCode(comIndex:comIndex+19)));
                % Maximal integration power
                maxPower = max(maxPower,comPower);
            end % Search different NH code combiniations
            fineResult(fineBinIndex) = maxPower;
        end % for numOfFineBins
        
        %--- Find the fine carrier freq. ----------------------------------
        [~, maxFinBin] = max(fineResult);
        acqResults.carrFreq(1,PRN) = fineFreqBins(maxFinBin);
        
        %signal found, if IF =0 just change to 1 Hz to allow processing
        if(acqResults.carrFreq(1,PRN) == 0)
            acqResults.carrFreq(1,PRN) = 1;
        end
        if (strcmp(acqType,'Normal'))==1
            %--- Indicate PRN number of the detected signal -------------------
            fprintf('Acquired %02d with Doppler %d Hz, Code Phase %d Samples and Test statistics %d \n'...
            , PRN, acqResults.carrFreq(1,PRN), acqResults.codePhase(1,PRN), acqResults.peakMetric(1,PRN));
        end
        %% Downsampling recovery ====================================
        % Find acquisition results corresponding to orignal sampling freq
        if (exist('oldFreq', 'var') && settings.resamplingflag == 1)
            % Find code phase
            acqResults.codePhase(1,PRN) = floor((primaryCodePhase - 1)/ ...
                settings.samplingFreq * oldFreq)+1;
            
            % Doppler frequency
            if (settings.IF >= settings.samplingFreq/2)
                % In this condition, the FFT computed freq. is symmetric
                % with the true frequemcy about half of the sampling
                % frequency, so we have the following:
                IF_temp = settings.samplingFreq - settings.IF;
                doppler = IF_temp - acqResults.carrFreq(1,PRN);
            else
                doppler = acqResults.carrFreq(1,PRN) - settings.IF;
            end
            
            % Carrier freq. corresponding to orignal sampling freq
            acqResults.carrFreq(1,PRN) = doppler + oldIF;
        end
        %% In case of APT acquisition, store the second peak
        if settings.AptActive==1  
            %--- Store result secondary peak in case it overcomes the APT threshold-----------------------------------------------------
            %acqResults.peakMetric(2,PRN) = secondaryPeakSize/sigPower/settings.acqNonCohTime;
            acqResults.peakMetric(2,PRN) = secondaryPeakSize/thirdPeakSize;
            acqResults.carrFreq(2,PRN) = fineFreqBins(maxFinBin);
            acqResults.codePhase(2,PRN) =secondaryPeakCodePhase;%secondaryPeakCodePhase; 
            
            if (strcmp(acqType,'APT')==1)
                if (acqResults.peakMetric(2,PRN)>settings.AptThreshold)
                    disp(['--------------SPOOFING STATUS PRN', num2str(PRN), ': SPOOFING DETECTED---------------'])
                    disp(['     Primary peak:   magnitude =', num2str(round(acqResults.peakMetric(1,PRN),2)) ,' Code-phase =', num2str(acqResults.codePhase(1,PRN))])
                    disp(['     Secondary peak: magnitude =', num2str(round(acqResults.peakMetric(2,PRN),2)) ,' Code-phase =', num2str(acqResults.codePhase(2,PRN))])
                    %disp('--------------------------------------')
                else
                    disp(['--------------SPOOFING STATUS PRN ', num2str(PRN), ': No spoofing---------------'])
                    disp(['     Primary peak magnitude =', num2str(acqResults.peakMetric(1,PRN)) ,' Code-phase =', num2str(acqResults.codePhase(1,PRN))])
                    disp(['     Secondary peak magnitude =', num2str(acqResults.peakMetric(2,PRN)) ,' Code-phase =', num2str(acqResults.codePhase(2,PRN))])
                    %disp('--------------------------------------')
                end
            end
            
            
        end
    else
        %--- No signal with this PRN --------------------------------------
        %fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%=== Acquisition is over ==================================================
fprintf(' \n');
