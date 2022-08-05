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