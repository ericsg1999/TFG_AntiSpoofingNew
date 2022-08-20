function [FirstSubFrameIndex]=detectFirstSubFrameIndex(bits, D30Star)

    %--- "Cut" one sub-frame's bits ---------------------------------------
    subframe = bits(1 : 300);

    %--- Correct polarity of the data bits in all 10 words ----------------
    for j = 1:10
        [subframe(30*(j-1)+1 : 30*j)] = ...
            checkPhase(subframe(30*(j-1)+1 : 30*j), D30Star);
        
        D30Star = subframe(30*j);
    end

    %--- Decode the sub-frame id ------------------------------------------
    % For more details on sub-frame contents please refer to GPS IS.
    FirstSubFrameIndex = bin2dec(subframe(50:52));
    
end

