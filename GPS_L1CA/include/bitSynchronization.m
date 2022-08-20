function [navBitsBin] = bitSynchronization(navBitsSamples)

%Bit synchronization is the process of adapting the bit rate coming from the tracking module to the PVT module.
% Due to the noise and weak power of the signal, a correlator output value
% shall be "averaged" among other consecutive values in order to mitigate the noise effect. In particular, 20
% successive correlator output values are "averaged" when deciding a PVT value. As a result
% of replacing 20 samples by a single one, the sample rate is reduced from 1KSPS (1 sample
% each ms) to 50SPS. 

%--- Group every 20 vales of bits into columns ------------------------
navBitsSamples = reshape(navBitsSamples, ...
    20, floor((size(navBitsSamples, 1) / 20)));

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

end

