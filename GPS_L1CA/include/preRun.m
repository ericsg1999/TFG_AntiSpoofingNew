% function [channel] = preRun(acqResults, settings)
% %Function initializes tracking channels from acquisition data. The acquired
% %signals are sorted according to the signal strength. This function can be
% %modified to use other satellite selection algorithms or to introduce
% %acquired signal properties offsets for testing purposes.
% %
% %[channel] = preRun(acqResults, settings)
% %
% %   Inputs:
% %       acqResults  - results from acquisition.
% %       settings    - receiver settings
% %
% %   Outputs:
% %       channel     - structure contains information for each channel (like
% %                   properties of the tracked signal, channel status etc.). 
% 
% %--------------------------------------------------------------------------
% %                           SoftGNSS v3.0
% % 
% % Copyright (C) Darius Plausinaitis
% % Written by Darius Plausinaitis
% % Based on Peter Rinder and Nicolaj Bertelsen
% %--------------------------------------------------------------------------
% %This program is free software; you can redistribute it and/or
% %modify it under the terms of the GNU General Public License
% %as published by the Free Software Foundation; either version 2
% %of the License, or (at your option) any later version.
% %
% %This program is distributed in the hope that it will be useful,
% %but WITHOUT ANY WARRANTY; without even the implied warranty of
% %MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% %GNU General Public License for more details.
% %
% %You should have received a copy of the GNU General Public License
% %along with this program; if not, write to the Free Software
% %Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
% %USA.
% %--------------------------------------------------------------------------
% 
% %CVS record:
% %$Id: preRun.m,v 1.8.2.20 2006/08/14 11:38:22 dpl Exp $
% 
% %% Initialize all channels ================================================
% channel                 = [];   % Clear, create the structure
% 
% channel.PRN             = 0;    % PRN number of the tracked satellite
% channel.acquiredFreq    = 0;    % Used as the center frequency of the NCO
% channel.codePhase       = 0;    % Position of the C/A  start
% 
% channel.status          = '-';  % Mode/status of the tracking channel
%                                 % "-" - "off" - no signal to track
%                                 % "T" - Tracking state
% 
% %--- Copy initial data to all channels ------------------------------------
% channel = repmat(channel, 1, settings.numberOfChannels);
% 
% %% Copy acquisition results ===============================================
% 
% %--- Sort peaks to find strongest signals, keep the peak index information
% [junk, PRNindexes]          = sort(acqResults.peakMetric, 2, 'descend');
% 
% %--- Load information about each satellite --------------------------------
% % Maximum number of initialized channels is number of detected signals, but
% % not more as the number of channels specified in the settings.
% for ii = 1:min([settings.numberOfChannels, sum(acqResults.carrFreq ~= 0)])
%     channel(ii).PRN          = PRNindexes(ii);
%     channel(ii).acquiredFreq = acqResults.carrFreq(PRNindexes(ii));
%     channel(ii).codePhase    = acqResults.codePhase(PRNindexes(ii));
%     
%     % Set tracking into mode (there can be more modes if needed e.g. pull-in)
%     channel(ii).status       = 'T';
% end
function [channel] = preRun(acqResults, settings)
%Function initializes tracking channels from acquisition data. The acquired
%signals are sorted according to the signal strength. This function can be
%modified to use other satellite selection algorithms or to introduce
%acquired signal properties offsets for testing purposes.
%
%[channel] = preRun(acqResults, settings)
%
%   Inputs:
%       acqResults  - results from acquisition.
%       settings    - receiver settings
%
%   Outputs:
%       channel     - structure contains information for each channel (like
%                   properties of the tracked signal, channel status etc.). 

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
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
%$Id: preRun.m,v 1.8.2.20 2006/08/14 11:38:22 dpl Exp $

%% Initialize all channels ================================================
channel                 = [];   % Clear, create the structure

channel.PRN             = 0;    % PRN number of the tracked satellite
channel.acquiredFreq    = 0;    % Used as the center frequency of the NCO
channel.codePhase       = 0;    % Position of the C/A  start

channel.status          = '-';  % Mode/status of the tracking channel
                                % "-" - "off" - no signal to track
                                % "T" - Tracking state

%--- Copy initial data to all channels ------------------------------------
% channel will have numSatToProcess*numberChannelsPerSat rows. Let's
% find numSatToProcess and numberChannelsPerSat

% How many satellites will be processed? It will be processed the detected satellites. However, if the number of detected satellites exceeds the 
% maximum number of satellites to process (settings.maxNumSatToProcess), some detected satellites will be dismissed. 
numSatToProcess=min([settings.maxNumSatToProcess, sum(acqResults.carrFreq(1,:) ~= 0)]);

%if APT spoofing detection is active, it will be allocated more than one
%channel (specified in settings.numberChannelsPerSat) per satellite. In
%case APT spoofing detection is not active, only one channel per satellite
%is enough.
numberChannelsPerSat=0;
if settings.AptActive==1
    numberChannelsPerSat=settings.AptNumberChannelsPerSat;
else
    numberChannelsPerSat=1;
end                               

channel = repmat(channel, 1, numSatToProcess*numberChannelsPerSat);

%% Copy acquisition results ===============================================

%--- Sort peaks to find strongest signals, keep the peak index information
% la funcio sort agafa la primera fila (peak1) de la matriu acqResults.peakMetric i la ordena 
% unica fila (per aixo el 2, sense el 2 ordena les columnes) en ordre
% descendent. junk es un vector igual que acqResults.peakMetric amb els
% valors ordenats (no ens interessa) i PRN indexes es un vector amb els
% index dels satel.lits amb el peakMetric de gran a petit
[junk, PRNindexes]          = sort(acqResults.peakMetric(1,:), 2, 'descend');


 %--- Load information about each satellite --------------------------------
channel_count=1;
for i = 1:numSatToProcess
    for j =1:numberChannelsPerSat
        channel(channel_count).PRN          = PRNindexes(i);
        if (j==1)
            channel(channel_count).acquiredFreq = acqResults.carrFreq(1,PRNindexes(i));
            channel(channel_count).codePhase    = acqResults.codePhase(1,PRNindexes(i));
            % Set tracking into mode (there can be more modes if needed e.g. pull-in)
            channel(channel_count).status       = 'T';  
        else
            channel(channel_count).acquiredFreq = acqResults.carrFreq(2,PRNindexes(i));
            channel(channel_count).codePhase    = acqResults.codePhase(2,PRNindexes(i));
            channel(channel_count).status       = 'APT';  
        end
        channel_count=channel_count+1;
    end
end