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

%%
I_P_InputBits=1:1:1000000;

subFrameStart=21;

navBitsSamples = I_P_InputBits(subFrameStart - 20 : ...
    subFrameStart + (1500 * 20) -1)';


prova=I_P_InputBits(subFrameStart+1500*20-20 : subFrameStart + 2*(1500 * 20) -1)';
prova2=I_P_InputBits(subFrameStart+2*1500*20-20 : subFrameStart + 3*(1500 * 20) -1)';

%%
a=[]
a=[a; [1 2 3 4 5]]
a=[a; [1 2 3 4 5]]
 %% Look for correlation peaks for coarse acquisition ============

%  results=[1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 ; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2;...
%           1 2 3 60 90 30 2 3 8 2 3 4 20 50 21 3 4 1 2 3 4 2 5 3 2 1 8 3 4 5 4 1 3 2; ...
%           1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2]
%  results=[1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 ; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2;...
%           90 30 1 2 3 2 3 8 2 3 4 20 50 21 3 4 60 2 3 4 2 5 3 2 1 8 3 4 5 4 1 3 2 1; ...
%           1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2]
%  
 results=[1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 ; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2;...
          61 2 1 2 8 2 1 2 3 30 50 31 1 1 2 60 90 1 2 3 2 1 88 2 1 2 3 4 2 1 2 3 1 3; ...
          1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2]
 
%     % Find the correlation peak and the carrier frequency
%     [~, acqCoarseBin] = max(max(results, [], 2));
%     % Find code phase of the same correlation peak
%     [peakSize, codePhase] = max(max(results));
%      
    
    
    %--------------- FIRST PEAK -------------------------------------------
    %--- Find the correlation peak and the carrier frequency --------------
    %M = max(A,[],dim) returns the maximum element along dimension dim. For example, if A is a matrix, then max(A,[],2) is a column vector containing the maximum value of each row.
    [primaryPeakSize frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [primaryPeakSize primaryCodePhase] = max(max(results));
    
    codephase_inferiorLimit=1;
    codephase_superiorLimit=size(results,2);
    samplesPerCode=17;
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
    samplesPerCodeChip=2;
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
    
    

%% NEW
I_P_InputBits=1:1000000;
subFrameStart=100;
InputBitsSamples=size(I_P_InputBits,2); %number of samples of the incoming bits
fiveSubframesSamples=1500*20; %number of bits of each frame (or 5 subframes)

numFiveSubframesInInputBits=floor((InputBitsSamples-subFrameStart)/fiveSubframesSamples);%number of frames inside input bits. It starts to count when it finds the subFrameStart
it=0;
eph = eph_structure_init();
%The next loop goes through all I_P_InputBits taking groups of 5 subframes
%(including the previous last bit) and decodes its ephemeris
for i=1:fiveSubframesSamples:numFiveSubframesInInputBits*fiveSubframesSamples %i=1:fiveSubframesSamples:InputBitsSamples
    
    navBitsSamples=I_P_InputBits(subFrameStart+it*(1500*20)-20 : subFrameStart + (it+1)*(1500 * 20) -1)';
    it=it+1;
    
    %--- Group every 20 vales of bits into columns ------------------------
    navBitsSamples = reshape(navBitsSamples, ...
        20, (size(navBitsSamples, 1) / 20));

    %--- Sum all samples in the bits to get the best estimate -------------
    navBits = sum(navBitsSamples);

    %--- Now threshold and make 1 and 0 -----------------------------------
    % The expression (navBits > 0) returns an array with elements set to 1
    % if the condition is met and set to 0 if it is not met.
    navBits = (navBits > 0);

    %--- Convert from decimal to binary -----------------------------------
    % The function ephemeris expects input in binary form. In Matlab it is
    % a string array containing only "0" and "1" characters.
    navBitsBin = dec2bin(navBits);
    
    
    %if NAVI --> ephemeris_new
    % if I only want to obtain position--> OLD
    %=== Decode ephemerides and TOW of the first sub-frame ================

%     [eph] = ephemeris_new(navBitsBin(2:1501)', navBitsBin(1),settings,eph);
    navBitsBin(2:1501)';
      
end
[naviTowDetectionResult] = naviTowSpoofingDetection(eph)

%% NEW 2.0
I_P_InputBits=1:1000000;
subFrameStart=100;

eph = eph_structure_init();
%The next loop goes through all I_P_InputBits taking groups of 5 subframes
%(including the previous last bit) and decodes its ephemeris

%Let's process the first frame, which may be not include all 5 subframes

navBitsSamplesFirstFiveSubFrames=I_P_InputBits(subFrameStart-20 : subFrameStart + (1500 * 20) -1)'; %The first frame could be incomplete. Let's take the first 5 subframes which can be of different frames

navBitsBinFirstFiveSubFrames = bitSynchronization(navBitsSamplesFirstFiveSubFrames); %In bits (1 bit = 20 samples)
    
%FirstSubFrameIndex=detectFirstSubFrameIndex(navBitsBinFirstFiveSubFrames(2:301)', navBitsBinFirstFiveSubFrames(1)); %ID of the first subframe received
FirstSubFrameIndex=3;     
numSubFramesFirstFrame=5-FirstSubFrameIndex+1; %How many subframes contains the first received frame

lastBitOfFirstFrame=numSubFramesFirstFrame*300+1;

% [eph] = ephemeris_new(navBitsBinFirstFiveSubFrames(2:lastBitOfFirstFrame)', navBitsBinFirstFiveSubFrames(1),settings,eph);
a=navBitsBinFirstFiveSubFrames(2:lastBitOfFirstFrame)';

fiveSubframesSamples=1500*20; %number of samples of each frame (or 5 subframes)

InputBitsSamples=size(I_P_InputBits,2); %number of samples of the incoming bits

numFiveSubframesInInputBits=floor((InputBitsSamples-subFrameStart-lastBitOfFirstFrame)/fiveSubframesSamples);%number of frames inside input bits. It starts to count when it finds the subFrameStart

it=0;

%Now that we have already decoded the first and problematic frame (as it
%may be incomplete), let's decode all frames 
for i=1:fiveSubframesSamples:numFiveSubframesInInputBits*fiveSubframesSamples %i=1:fiveSubframesSamples:InputBitsSamples
    
    %navBitsSamples=I_P_InputBits(subFrameStart+it*(1500*20)-20 : subFrameStart + (it+1)*(1500 * 20) -1)';
    
%     firstCompleteFrameSample=subFrameStart+numSubFramesFirstFrame*300*20;
    
%     firstCompleteFrameSample=subFrameStart+lastBitOfFirstFrame*20+it*1500;
    firstCompleteFrameSample=subFrameStart+lastBitOfFirstFrame*20;
    firstFrameSample=firstCompleteFrameSample+(it)*(1500*20)-20;
    %firstInputBit=subFrameStart+numSubFramesFirstFrame*300*20+it*(1500*20)-20;
    lastFrameSample=firstFrameSample +20+ (1500 * 20)-1;
    navBitsSamples=I_P_InputBits(firstFrameSample : lastFrameSample)';
    
    
    it=it+1;
    
    navBitsBin = bitSynchronization(navBitsSamples);
    
    
    
%     if i==1 %PROCESSING OF THE FIRST FRAME. Sólo le pasaremos los subframes del primer frame. Por ejemplo, si el primer subframe que recibimos es el 3, la primera vez solo decodificaremos les efemerides de los subframes 3,4 y 5
%          %when p
%         
%     else %PROCESSING DE LOS SUBFRAMES QUE NO SEAN EL PRIMERO. LE PASASMOS 5 SUBFRAMES
%         firstBitOfFrame=numSubFramesFirstFrame*300+1+1500*(i-1);
%         
%         lastBitOfFirstFrame=firstBitOfFrame+1500;
        
%         [eph] = ephemeris_new(navBitsBin(2:1501)', navBitsBin(1),settings,eph); %when p
    b=navBitsBin(2:1501)';
end
[naviTowDetectionResult] = naviTowSpoofingDetection(eph)
 
%%
a=zeros(1,10)
a(1)=2
a(2)=4

 %% Look for correlation peaks for coarse acquisition ============

%  results=[1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 ; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2;...
%           1 2 3 60 90 30 2 3 8 2 3 4 20 50 21 3 4 1 2 3 4 2 5 3 2 1 8 3 4 5 4 1 3 2; ...
%           1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2]
%  results=[1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 ; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2;...
%           90 30 1 2 3 2 3 8 2 3 4 20 50 21 3 4 60 2 3 4 2 5 3 2 1 8 3 4 5 4 1 3 2 1; ...
%           1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2]
%  
%  results=[1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 ; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2;...
%           61 2 1 2 8 2 1 2 3 30 50 31 1 1 2 60 90 1 2 3 2 1 88 2 1 2 3 4 2 1 2 3 1 3; ...
%           1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2]
%   results=[1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 ; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2;...
%           10 30 1 2 3 2 3 8 2 3 4 20 50 21 3 4 60 2 3 37 41 5 3 2 1 80 90 70 5 4 1 33 2 1; ...
%           1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2]
  results=[1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 ; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2;...
          10 30 1 2 3 2 3 8 2 3 4 20 150 21 3 4 60 2 3 37 41 5 3 2 1 8 0 7 5 4 1 33 2 1; ...
          1 2 3 1 1 2 3 4 1 4 5 1 2 4 1 2 1 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2; 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2 1 2 3 4 2 5 3 2 1 2 3 4 5 4 1 3 2]
%     % Find the correlation peak and the carrier frequency
%     [~, acqCoarseBin] = max(max(results, [], 2));
%     % Find code phase of the same correlation peak
%     [peakSize, codePhase] = max(max(results));
%      
    
    
    %--------------- FIRST PEAK -------------------------------------------
    %--- Find the correlation peak and the carrier frequency --------------
    %M = max(A,[],dim) returns the maximum element along dimension dim. For example, if A is a matrix, then max(A,[],2) is a column vector containing the maximum value of each row.
    [primaryPeakSize frequencyBinIndex] = max(max(results, [], 2));

    %--- Find code phase of the same correlation peak ---------------------
    [primaryPeakSize primaryCodePhase] = max(max(results));
    
    codephase_inferiorLimit=1;
    codephase_superiorLimit=size(results,2);
    samplesPerCode=17;
    
    ms=0;%will indicate at which ms will be the first peak
    %if primaryCodePhase is in the 1st ms, look for the second and third peak in
    %the 1st ms. If primaryCodePhase is in the 2nd ms, look for the
    %second and third peak in the 2nd ms
    if primaryCodePhase<=samplesPerCode %the first peak is in the 1st ms
        ms=1;
        codephase_inferiorLimit=1;
        codephase_superiorLimit=samplesPerCode;
    else %the first peak is in the 2nd ms
        ms=2;
        codephase_inferiorLimit=samplesPerCode+1;
        codephase_superiorLimit=size(results,2);
    end
   
    %--------------- SECOND PEAK ------------------------------------------
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
%     samplesPerCodeChip   = round(settings.samplingFreq / settings.codeFreqBasis);
    samplesPerCodeChip=2;
    excludeRangeIndex1 = primaryCodePhase - samplesPerCodeChip;
    excludeRangeIndex2 = primaryCodePhase + samplesPerCodeChip;
    
    caso=0;
    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < codephase_inferiorLimit+2
        caso=1;
        codePhaseRange = excludeRangeIndex2 : ...
                         (codephase_superiorLimit-1 + excludeRangeIndex1);
                     
    elseif excludeRangeIndex2 >= codephase_superiorLimit
        caso=2;
        codePhaseRange = (excludeRangeIndex2 - samplesPerCode) : ...
                         excludeRangeIndex1;
    else
        caso=3;
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

    
    if caso==1
        if excludeRangeIndex11<=excludeRangeIndex2
            codePhaseRange2=[excludeRangeIndex22:codephase_superiorLimit];
        elseif excludeRangeIndex22>=codephase_superiorLimit
            codePhaseRange2=[excludeRangeIndex2:excludeRangeIndex11];
        else
            codePhaseRange2=[excludeRangeIndex2:excludeRangeIndex11, excludeRangeIndex22:codephase_superiorLimit];
        end
            

        
    elseif caso==2
        if excludeRangeIndex11<=codephase_inferiorLimit
            codePhaseRange2=[excludeRangeIndex22:codephase_superiorLimit];
        elseif excludeRangeIndex22>=excludeRangeIndex1
            codePhaseRange2=[codephase_inferiorLimit:excludeRangeIndex11];
        else
            codePhaseRange2=[codephase_inferiorLimit:excludeRangeIndex11, excludeRangeIndex22:excludeRangeIndex1];
        end
        
        
    elseif caso==3
        if excludeRangeIndex11<=codephase_inferiorLimit
            codePhaseRange2=[excludeRangeIndex22:excludeRangeIndex1,excludeRangeIndex2:codephase_superiorLimit];
        elseif excludeRangeIndex22>=excludeRangeIndex1
            codePhaseRange2=[codephase_inferiorLimit:excludeRangeIndex1,excludeRangeIndex2:excludeRangeIndex11];
        else
            if excludeRangeIndex1<excludeRangeIndex11 %primary peak is at the left of the secondary peak
                codePhaseRange2=[codephase_inferiorLimit:excludeRangeIndex1, excludeRangeIndex2:excludeRangeIndex11, excludeRangeIndex22:codephase_superiorLimit];
            else %primary peak is at the right of the secondary peak
                codePhaseRange2=[codephase_inferiorLimit:excludeRangeIndex11, excludeRangeIndex22:excludeRangeIndex1, excludeRangeIndex2:codephase_superiorLimit];
            end
            
        end
  
    end
        
%     %--- Find the third highest correlation peak in the same freq. bin ---
    [thirdPeakSize junk] = max(results(frequencyBinIndex, codePhaseRange2));
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
%     %--- Correct C/A code phase exclude range if the range includes array
%     %boundaries
%     if excludeRangeIndex1 < 2
%         
%         codePhaseRange2 = [excludeRangeIndex2 : excludeRangeIndex11...
%                           excludeRangeIndex22:(codephase_superiorLimit-1 + excludeRangeIndex1)];
% 
%     elseif excludeRangeIndex2 >= samplesPerCode
%         codePhaseRange2 = [(excludeRangeIndex2 - samplesPerCode):excludeRangeIndex11 ...
%                          excludeRangeIndex22:excludeRangeIndex1] ;
%     else
%         if excludeRangeIndex1<excludeRangeIndex11
%             codePhaseRange2 = [codephase_inferiorLimit:excludeRangeIndex1  excludeRangeIndex2:excludeRangeIndex11 excludeRangeIndex22:codephase_superiorLimit];
%         else
%             codePhaseRange2 = [codephase_inferiorLimit:excludeRangeIndex11  excludeRangeIndex22:excludeRangeIndex1 excludeRangeIndex2:codephase_superiorLimit];
%         end
%     end
%     
%     %--- Find the third highest correlation peak in the same freq. bin ---
%     [thirdPeakSize junk] = max(results(frequencyBinIndex, codePhaseRange2));
%     
    
%% 
a= [1 2 3 4 5 1; 3 4 5 3 2 1; 3 5 3 29 191 1]
b=[1 2 3 4 50 6 7 80 9 10]
c=[b(1:4) b(7:9)]

max(b(c))
find(b==max(b(1,c)))

%%
a= [-1 -1 -1 1 -1 -1 -1 1 1 1]
negative=find(a==-1)
positive=find(a==1)
size(negative,2)
size(positive,2)
size(a,2)