function [naviTowDetectionResult] = naviTowSpoofingDetection(eph)

naviTowDetectionResult=0; % 0 if no spoofing, otherwise indicates the TOW where the possible spoofing attack is occuring
framesToProcess=size(eph.TOW,1);
subframesPerFrame=5; %each frame contains 5 subframes

firstFrameTows=eph.TOW(1,:)
frameMinimumTow=min(firstFrameTows(firstFrameTows>0)); 
%subframeIndexStart=find(eph.TOW(1,:)==frameMinimumTow); %finds which is the first subframe (out of the 5 subframes) that time-wise was received first

% previousTOW=0;
% actualTOW=eph.TOW(1,subframeIndexStart);
previousTOW=eph.TOW(1,1);
actualTOW=eph.TOW(1,2);

for i=1:1:framesToProcess
    
    frameTransition=0;
    
    for j=1:1:subframesPerFrame
        if((j==subframesPerFrame)&&(i~=framesToProcess))
            frameTransition=1;
            previousTOW=eph.TOW(i,j);
            actualTOW=eph.TOW(i+1,1);
        elseif((j==subframesPerFrame)&&(i==framesToProcess))
            break
        else
            previousTOW=eph.TOW(i,j);
            actualTOW=eph.TOW(i,j+1);
        end
       
        if ((previousTOW~=0)&&(actualTOW~=0))
            if (actualTOW-previousTOW==6)
                if frameTransition==0
                    disp(['NAVI: TOW consistency in frame number: ', num2str(i),'--> subframes number: ', num2str(j), ' and ', num2str(j+1), ' no SPOOFING'])
                else
                    disp(['NAVI: TOW consistency in frames number: ', num2str(i),' and ',num2str(i+1) ,' --> subframes number: ', num2str(j), ' and ', num2str(1), ' no SPOOFING'])
                end
                
                
            else
                if frameTransition==0
                    disp('############### NAVI ALERT #########################')
                    disp(['NAVI: TOW inconsistency in frame number: ', num2str(i),'--> subframes number: ', num2str(j), ' and ', num2str(j+1)])
                else
                    disp(['NAVI: TOW inconsistency in frames number: ', num2str(i),' and ',num2str(i+1) ,' --> subframes number: ', num2str(j), ' and ', num2str(1)])
                end
                
            end
        else
            %Do nothing, don't check TOW consistency
        end
        
    end
end
% for i=1:1:framesToProcess
%     
%     actualTOW=eph.TOW(i,subframeIndexStart);
%     
%     for j=subframeIndexStart:1:(subframeIndexStart+subframesPerFrame-1)
%         
%         if j<=subframesPerFrame
%             previousSubframe=j-1;
%             actualSubframe=j;
%         elseif j==subframesPerFrame+1
%             previousSubframe=j-1;
%             actualSubframe=j-subframesPerFrame;
%         else
%             previousSubframe=j-subframesPerFrame;
%             actualSubframe=j-subframesPerFrame+1;
%         end
%         
%         if ((i==1 && j==subframeIndexStart)==0) 
%             if (actualTOW-previousTOW==6)
%                 if j==subframeIndexStart
%                     disp(['NAVI: TOW consistency in frame numbers: ', num2str(i-1),' and ', num2str(i), '--> subframes number: ', num2str(previousSubframe), ' and ', num2str(actualSubframe), ' no SPOOFING'])
%                 else
%                     disp(['NAVI: TOW consistency in frame number: ', num2str(i), '--> subframes number: ', num2str(previousSubframe), ' and ', num2str(actualSubframe), ' no SPOOFING'])
%                 end
%                 
% %                 disp(['NAVI: TOW consistency in frame number: ', num2str(i), '--> subframes number: ', num2str(j-1), ' and ', num2str(j), ' no SPOOFING'])
%             else
%                 disp('############### NAVI ALERT #########################')
%                 disp(['NAVI: TOW inconsistency in frame number: ', num2str(i), '--> subframes number: ', num2str(previousSubframe), ' and ', num2str(actualSubframe), 'possible SPOOFING'])
%             end
%         end
% 
%         if j<subframesPerFrame
%             previousTOW=eph.TOW(i,j);
%             actualTOW=eph.TOW(i,j+1);
%         elseif j==subframesPerFrame
%             previousTOW=eph.TOW(i,j);
%             actualTOW=eph.TOW(i,j-subframesPerFrame+1);
%         else
%             previousTOW=eph.TOW(i,j-subframesPerFrame);
%             actualTOW=eph.TOW(i,j-subframesPerFrame+1);
%             
%         end
%     end
% end

end

