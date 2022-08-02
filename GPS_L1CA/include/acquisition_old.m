function acqResults = acquisition_old(longSignal, settings, acqType,SatellitePresentList)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(longSignal, settings)
%
%   Inputs:
%               longSignal    - 11 ms of raw signal from the front-end 
%               settings      - Receiver settings. Provides information about
%                               sampling and intermediate frequencies and other
%                               parameters including the list of the satellites to
%                               be acquired.
%               acqType      - Describes if it is the 'NORMAL' acquisition (looking
%                               for present satellites) or 'APT' acquisition (looking for secondary peaks in only the present satellites alerting for spoofing).
%                               For the 'APT' acquisition, previously it shall be
%                               performed a 'NORMAL' acquisition that finds the
%                               present satellites, as 'APT' acquisition will only
%                               generate the grid for these present satellites (efficiency
%                               purposes). 
%       present_sat_PRN_list  - Vector containing the PRN of the present
%                              satellites. Recall that the acqResults
%                              (output of this function) has an attribute
%                              that is a list of present satellites' PRN
%   Outputs:
%               acqResults    - Function saves code phases and frequencies of the 
%                               detected signals in the "acqResults" structure. The
%                               field "carrFreq" is set to 0 if the signal is not
%                               detected for the given PRN number. 
 
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis and Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on Peter Rinder and Nicolaj Bertelsen
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

%% Initialization =========================================================
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
                        (settings.codeFreqBasis / settings.codeLength));

% Takes the first 2 msec  and store them in different vectors and one with zero DC
signal1 = longSignal(1 : samplesPerCode);
signal2 = longSignal(samplesPerCode+1 : 2*samplesPerCode);

 
signal0DC = longSignal;

% signal0DC = longSignal - mean(longSignal);

% Find sampling period
ts = 1 / settings.samplingFreq;

% Find phase points of the local carrier wave 
phasePoints = (0 : (samplesPerCode-1)) * 2 * pi * ts;

% Number of the frequency bins for the given acquisition band (500Hz steps)
numberOfFrqBins = round((settings.acqSearchBand)*1000/settings.acqSearchStep) + 1;

% Generate all C/A codes and sample them according to the sampling freq.
caCodesTable = makeCaTable_old(settings);

%--- Input signal power for GLRT statistic calculation --------------------
sigPower = sqrt(var(longSignal(1:samplesPerCode)) * samplesPerCode);

%--- Initialize arrays to speed up the code -------------------------------
% Search results of all frequency bins and code shifts (for one satellite)
results     = zeros(numberOfFrqBins, samplesPerCode);

% Carrier frequencies of the frequency bins
frqBins     = zeros(1, numberOfFrqBins);

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


%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);
acqResults.SatellitePresentList= [];
fprintf('(');

% Depending on the acqType, the acquisition will be performed among all the
% possible satellites(Normal acquisition) or only the present/visible ones
% (APT)
list_sat_to_acquire=zeros(1,32);
switch acqType
    case 'Normal'
        list_sat_to_acquire=settings.acqSatelliteList;
    case 'APT'
        list_sat_to_acquire=SatellitePresentList;
    otherwise 
        fprintf('ERROR: Incorrect Acquisition type (acqType) input parameter of the function acquisition.m');
end
    
% for PRN = settings.acqSatelliteList
for PRN = list_sat_to_acquire
%% Correlate signals ======================================================   
    %--- Perform DFT of C/A code ------------------------------------------
    caCodeFreqDom = conj(fft(caCodesTable(PRN, :)));

    %--- Make the correlation for whole frequency band (for all freq. bins)
    for frqBinIndex = 1:numberOfFrqBins

        %--- Generate carrier wave frequency grid (0.5kHz step) -----------
        frqBins(frqBinIndex) = settings.IF - ...
                               (settings.acqSearchBand/2) * 1000 + ...
                               settings.acqSearchStep * (frqBinIndex - 1);

        %
        switch settings.fileType
            case 1
                %--- Generate local sine and cosine -------------------------------
                sinCarr = sin(frqBins(frqBinIndex) * phasePoints);
                cosCarr = cos(frqBins(frqBinIndex) * phasePoints);

                %--- "Remove carrier" from the signal -----------------------------
                I1      = sinCarr .* signal1;
                Q1      = cosCarr .* signal1;
                I2      = sinCarr .* signal2;
                Q2      = cosCarr .* signal2;

                %--- Convert the baseband signal to frequency domain --------------
                IQfreqDom1 = fft(I1 + j*Q1);
                IQfreqDom2 = fft(I2 + j*Q2); 
            
            case 2
                %--- Generate local exponential -------------------------------
                expCarr=exp(-1j*frqBins(frqBinIndex) * phasePoints); %e^(-j*2pi*fd_estimation*t)
                %% INTEGRATION
                        %--- Do correlation -----------------------------------------------
                for nonCohIndex = 1: settings.acqNonCohTime
                    % Take 2ms vectors of input data to do correlation
                    signal = longSignal((nonCohIndex - 1) * samplesPerCode + ...
                        1 : (nonCohIndex) * samplesPerCode);
                    % "Remove carrier" from the signal
                    I      = real(expCarr .* signal);
                    Q      = imag(expCarr .* signal);
                    % Convert the baseband signal to frequency domain
                    IQfreqDom = fft(I + 1i*Q);
                    % Multiplication in the frequency domain (correlation in
                    % time domain)
                    convCodeIQ = IQfreqDom .* caCodeFreqDom;
                    % Perform inverse DFT and store correlation results
                    cohRresult = abs(ifft(convCodeIQ));
                    % Non-coherent integration
                    results(frqBinIndex, :) = results(frqBinIndex, :) + cohRresult;
                end % nonCohIndex = 1: settings.acqNonCohTime
                
%                 %% NO INTEGRATION
% %                 %--- "Remove carrier" from the signal -----------------------------
%                 baseband_samples_complex1=signal1.*expCarr;
%                 baseband_samples_complex2=signal2.*expCarr;
% %                 %--- Convert the baseband signal to frequency domain --------------
%                 IQfreqDom1 = fft(baseband_samples_complex1);
%                 IQfreqDom2 = fft(baseband_samples_complex2);
                
        end
        

        %--- Multiplication in the frequency domain (correlation in time
        %domain)
        convCodeIQ = IQfreqDom .* caCodeFreqDom;
        %convCodeIQ2 = IQfreqDom2 .* caCodeFreqDom;

        %--- Perform inverse DFT and store correlation results ------------
        acqRes1 = abs(ifft(convCodeIQ)) .^ 2;
        %acqRes2 = abs(ifft(convCodeIQ2)) .^ 2;
        acqRes2=0;
        %--- Check which msec had the greater power and save that, will
        %"blend" 1st and 2nd msec but will correct data bit issues
        if (max(acqRes1) > max(acqRes2))
            results(frqBinIndex, :) = acqRes1;
        else
            results(frqBinIndex, :) = acqRes2;
        end
        
    end % frqBinIndex = 1:numberOfFrqBins
 %% plot the grid of the acquired acquisition search grids
 
    if settings.acqInitialPlots==1 && (strcmp(acqType,'Normal')==1)
        Td=0:1:(samplesPerCode-1);
        figure;
        mesh(Td,frqBins,results);
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
        
        mesh(Td,frqBins,results);
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

%% Look for correlation peaks in the results ==============================
    % Find the highest peak and compare it to the second highest peak
    % The second peak is chosen not closer than 1 chip to the highest peak
    
    %--------------- FIRST PEAK -------------------------------------------
    %--- Find the correlation peak and the carrier frequency --------------
    %M = max(A,[],dim) returns the maximum element along dimension dim. For example, if A is a matrix, then max(A,[],2) is a column vector containing the maximum value of each row.
    [primaryPeakSize frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [primaryPeakSize primaryCodePhase] = max(max(results));
    
    
    %--------------- SECOND PEAK ------------------------------------------
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    excludeRangeIndex1 = primaryCodePhase - samplesPerCodeChip;
    excludeRangeIndex2 = primaryCodePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange = excludeRangeIndex2 : ...
                         (samplesPerCode + excludeRangeIndex1-1);
                         
    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
                          excludeRangeIndex2 : samplesPerCode];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    [secondaryPeakSize secondaryPeakCodePhase] = max(results(frequencyBinIndex, codePhaseRange));
    
    
    %--------------- THIRD PEAK -------------------------------------------
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    excludeRangeIndex11 = secondaryPeakCodePhase - samplesPerCodeChip;
    excludeRangeIndex22 = secondaryPeakCodePhase + samplesPerCodeChip;

    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 2
        codePhaseRange2 = [excludeRangeIndex2 : excludeRangeIndex11...
                          excludeRangeIndex22:(samplesPerCode + excludeRangeIndex1-1)];

    elseif excludeRangeIndex2 >= samplesPerCode
        codePhaseRange2 = [(excludeRangeIndex2 - samplesPerCode):excludeRangeIndex11 ...
                         excludeRangeIndex22:excludeRangeIndex1] ;
    else
        if excludeRangeIndex1<excludeRangeIndex11
            codePhaseRange2 = [1:excludeRangeIndex1  excludeRangeIndex2:excludeRangeIndex11 excludeRangeIndex22:samplesPerCode];
        else
            codePhaseRange2 = [1:excludeRangeIndex11  excludeRangeIndex22:excludeRangeIndex1 excludeRangeIndex2:samplesPerCode];
        end
    end

    %--- Find the third highest correlation peak in the same freq. bin ---
    [thirdPeakSize junk] = max(results(frequencyBinIndex, codePhaseRange2));
    
    %--- Store result -----------------------------------------------------
    acqResults.peakMetric(1,PRN) = primaryPeakSize/thirdPeakSize;            
    % If the result is above threshold, then there is a signal ...
%     if (primaryPeakSize/thirdPeakSize) >= settings.acqThreshold
%        acqResults.SatellitePresentList=[acqResults.SatellitePresentList PRN];
%        switch settings.acqFineFreqSearch
%            case 0
%                %--- Save properties of the detected satellite signal -------------
%                 acqResults.carrFreq(1,PRN)  = frqBins(frequencyBinIndex);
%                 acqResults.codePhase(1,PRN) = primaryCodePhase;
%            case 1
%                if frqBins(frequencyBinIndex)>0
%                    %% Fine resolution frequency search =======================================
% 
%                     %--- Indicate PRN number of the detected signal -------------------
%                     fprintf('%02d ', PRN);
% 
%                     %--- Generate 10msec long C/A codes sequence for given PRN --------
%                     caCode = generateCAcode(PRN);
% 
%                     codeValueIndex = floor((ts * (1:10*samplesPerCode)) / ...
%                                            (1/settings.codeFreqBasis));
% 
%                     longCaCode = caCode((rem(codeValueIndex, 1023) + 1));
% 
%                     %--- Remove C/A code modulation from the original signal ----------
%                     % (Using detected C/A code phase)
%                     xCarrier = ...
%                         signal0DC(primaryCodePhase:(primaryCodePhase + 10*samplesPerCode-1)) ...
%                         .* longCaCode;
% 
%                     %--- Find the next highest power of two and increase by 8x --------
%                     fftNumPts = 8*(2^(nextpow2(length(xCarrier))));
% 
%                     %--- Compute the magnitude of the FFT, find maximum and the
%                     %associated carrier frequency 
%                     fftxc = abs(fft(xCarrier, fftNumPts)); 
%                     %fd_axis=-7000+settings.FI:14000/length(fftxc):7000-14000/length(fftxc)+settings.FI;
%                     %plot(fftxc);
%                     uniqFftPts = ceil((fftNumPts + 1) / 2);
%                     [fftMax, fftMaxIndex] = max(fftxc(5 : uniqFftPts-5));
% 
%                     fftFreqBins = (0 : uniqFftPts-1) * settings.samplingFreq/fftNumPts;
% 
%                     %--- Save properties of the detected satellite signal -------------
%                     acqResults.carrFreq(1,PRN)  = fftFreqBins(fftMaxIndex);
%                     acqResults.codePhase(1,PRN) = primaryCodePhase;
%                else
%                     acqResults.carrFreq(1,PRN)  = -1;
%                     acqResults.codePhase(1,PRN) = -1;
%                end
% 
%            otherwise
%                disp('initSettings acqFineFreqSearch parameter holds an invalid value');
%% Fine carrier frequency search =================================================
    %--- Do fine acquisition -----------------------------------
    if acqResults.peakMetric(1,PRN) > settings.acqThreshold
        acqResults.SatellitePresentList=[acqResults.SatellitePresentList PRN];
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);
        
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
            fineFreqBins(fineBinIndex) = frqBins(frequencyBinIndex) + ...
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
        % Save code phase acquisition result
        acqResults.codePhase(1,PRN) = primaryCodePhase;
        %signal found, if IF =0 just change to 1 Hz to allow processing
        if(acqResults.carrFreq(1,PRN) == 0)
            acqResults.carrFreq(1,PRN) = 1;
        end
       

        
        %--- Indicate PRN number of the detected signal -------------------
        fprintf('Acquired %02d with Doppler %d Hz, Code Phase %d Samples and Test statistics %d \n'...
        , PRN, acqResults.carrFreq(1,PRN), acqResults.codePhase(1,PRN), acqResults.peakMetric(1,PRN));
    
    
    %% In case of APT acquisition, search of the second peak
        if settings.AptActive==1  
            %--- Store result secondary peak in case it overcomes the APT threshold-----------------------------------------------------
            acqResults.peakMetric(2,PRN) = secondaryPeakSize/thirdPeakSize;
            acqResults.carrFreq(2,PRN) = fineFreqBins(maxFinBin);
            acqResults.codePhase(2,PRN) =secondaryPeakCodePhase; 
                
        end

    else
        %--- No signal with this PRN --------------------------------------
        %fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList
%update the list of the visible satellites in the settings structure

%=== Acquisition is over ==================================================
fprintf(')\n');
end
