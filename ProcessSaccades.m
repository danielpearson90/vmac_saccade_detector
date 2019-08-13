function ProcessSaccades(minSub, maxSub, ROIs, varargin)
% ProcessSaccades: I-VT Algorithm for determining saccade direction and
% latency from Tobii eyetracker data (updated to work with Tobii Pro SDK)
%
% ProcessSaccades(subjectlist, ROIs, [fixationCoords], [discardAnticipatorySaccades], [discardOutsideFixationSaccades], [graphVersion])
% subjectlist is a vector of participant numbers, e.g., [1:5 8 9:15].
%
% ROIs is a list of the coordinates of stimuli of interest (in PTB pixel
% format), e.g., [0 0 500 500; 500 500 1000 1000]. If ROIs is set to 1, assumes saccades will be separated according to
% whether they are directed left or right relative to the fixation point.
% If ROIs set to 2, assumes saccades will be separated according to whether
% they are directed up or down relative to the fixation point.
%
% fixationCoords is an optional variable that codes for the participant's
% starting fixation point on each trial. If left out, algorith assumes the
% centre of the screen.
%
% discardAnticipatorySaccades is optional. If set to 0, anticipatory saccades will not be
% discarded. Default value is to discard anticipatory saccades.
%
% discardOutsideFixationSaccades is optional. If set to 0, trials with no
% valid gaze samples recorded within 100 px of the fixation point in first
% 80 ms will not be discarded. Default value is to discard.

% Add folders to path
addpath(genpath('sample_data'))
addpath(genpath('utility_functions'))

% Screen details
res = [1920 1080];

if nargin<1
    error('missing inputs')
elseif nargin > 7
    error('too many inputs')
end

nVarargs = length(varargin);
if nVarargs < 3
    varargin{3} = [];
end
if isempty(varargin{1})
    varargin{1} = [res(1)/2 res(2)/2];
end
if isempty(varargin{2})
    varargin{2} = 1;
end
if isempty(varargin{3})
    varargin{3} = 1;
end
if isempty(varargin{4})
	varargin{4} = false;
end
fixationCoords = varargin{1};
discardAnticipatory = varargin{2};
discardOutsideFixation = varargin{3};
graphVersion = varargin{4};



% Default Algorithm parameters -------------------------------------------------------------------
WindowLength = 20; % ms used in sampling window to determine angular velocity
vThreshold = 40; % velocity threshold (dva/s) for classification as a saccade

% if we do not include ROIs, then assume default VMAC locations
if isempty(ROIs)
    ROIs = [914.5, 294.5, 1005.5, 385.5;
          741.294919243112, 394.5, 832.294919243112, 485.5;
          741.294919243112, 594.5, 832.294919243112, 685.5;
          914.5, 694.5, 1005.5, 785.5;
          1087.70508075689, 594.5, 1178.70508075689, 685.5;
          1087.70508075689, 394.5, 1178.70508075689, 485.5];
end

fixGapDurThresh = 10; % a saccade needs to be longer than 10ms to be counted, otherwise it is likely a blink


%%

if discardAnticipatory == 1
    anticipationThreshold = 80; % saccades that are initiated <80 ms after trial begins will be discarded
else
    anticipationThreshold = 0;
end
%---------------------------------------------------------------------------------------------------


%% PARAMETERS TO BE UPDATED FOR EACH ANALYSIS -------------------------------------------------------
eyeDataFolder = 'sample_data/EyeData'; % path to folder containing raw gaze data
behavDataFolder = 'sample_data/BehavData'; % path to folder containing behavioural data for VMAC expt

eyeDataBase = '/GazeDataP'; % base filename for eye data
behavDataBase = '/VMAC_ET_xxx_P'; % base filename for behavioural data

maxPhases = [2]; % phases in the experiment, including the practice phase
maxTrials = [0, 60]; %number of trials to analyse for each phase (set first value to 0 to avoid analysing practice data)
maxSessions = 1; % number of sessions in the experiment

if graphVersion
	trialsToCheck = [1:480]; % input a vector of trial numbers here for graphs to be generated
end

%---------------------------------------------------------------------------------------------------
%% Make directory for processed data ---------------------------------------------------------------

mkdir(eyeDataFolder, 'processed_saccade_data');
mkdir(eyeDataFolder, 'raw_velocity_data');
mkdir(eyeDataFolder, 'all_eye_movements');

%%

%Calculate the centre point of each ROI
if numel(ROIs) > 1
    stimCentre = [(ROIs(:,3)+ROIs(:,1))/2 (ROIs(:,4)+ROIs(:,2))/2];
    withinAngle = 30; % max degrees of polar angle from stimulus for saccade to be classified as towards the stim.
elseif ROIs == 1
    stimCentre(1,:) = [1 res(2)/2]; %left location
    stimCentre(2,:) = [res(1) res(2)/2]; %right location
    withinAngle = 89; %set to whatever is appropriate, but be wary of 90 degs as classifier currently checks for >= this number
else
    stimCentre(1,:) = [res(1)/2 0]; %top location
    stimCentre(2,:) = [res(1)/2 res(2)]; %bottom location
    withinAngle = 89; %set to whatever is appropriate, but be wary of 90 degs as classifier currently checks for >= this number
end

%create vectors from the fixation point to the centre of the ROIs
stimVector = stimCentre - repmat(fixationCoords,size(stimCentre,1),1);

subStep = 0;
for sub = minSub : maxSub
    
    subStep = subStep+1; 
    
    for session = 1:maxSessions
        
        % This may need to be updated depending on filename structure
        GazeFileName = [eyeDataFolder, eyeDataBase, num2str(sub), 'S', num2str(session), '.mat']; 
        BehavFileName = [behavDataFolder, behavDataBase, num2str(sub), 'S', num2str(session), '.mat'];
        
        if exist(GazeFileName) == 2

	        load(GazeFileName, 'ALLGAZEDATA');
	        load(BehavFileName);
	        
	        for phase = 1:maxPhases(session)
	            
	            saccadeSessionData(session).saccadePhaseData(phase).summarySaccadeData = zeros(maxTrials(session, phase),size(stimVector,1)+9); %initiate summary data array
	            
	            clc; disp(['Processing Subject ',num2str(sub), ' Session ', num2str(session), ' Phase ', num2str(phase)])
	            
	            
	            
	            for t = 1:maxTrials(session, phase)
	                
	                clc; disp(['Sub ', num2str(sub), ' Session ', num2str(session), ' Phase ', num2str(phase), ' Trial ', num2str(t)]);
	                
	                %load data file for current trial
	                savedEGdata = ALLGAZEDATA.EGdataPhase(phase).EGdataTrial(t).data;
	                savedEGdata.system_time_stamp = double(savedEGdata.system_time_stamp);
	                
	                sampleInterval = ((savedEGdata.system_time_stamp(end,1)-savedEGdata.system_time_stamp(1,1))/1000)/size(savedEGdata.system_time_stamp,1); % determine the average interval between each successive eyetracker sample
	                sampleWindowLength = round((WindowLength/2)/sampleInterval); %determine how many samples for half of the saccade detection Window
	                
	                % Load some behav data for graphs later on
	                targetLoc = DATA.trialInfo(phase).trialData(t,5);
	                distractLoc = DATA.trialInfo(phase).trialData(t,6);


	                
	                % prepare essential EG data
	                EGerr = [mean(savedEGdata.left_gaze_point_validity) mean(savedEGdata.right_gaze_point_validity)]; % calc error on each eye
	                if EGerr(1) < EGerr(2)
	                	trialGazeData = RenameField(savedEGdata, {'left_gaze_origin_in_user_coordinate_system', 'left_gaze_origin_validity', 'left_gaze_point_in_user_coordinate_system', 'left_gaze_point_validity', 'left_gaze_point_on_display_area'},...
	                		{'gaze_origin_3d', 'gaze_origin_validity', 'gaze_point_3d', 'gaze_point_validity', 'gaze_point_2d'});

	                else
	                    trialGazeData = RenameField(savedEGdata, {'right_gaze_origin_in_user_coordinate_system', 'left_gaze_origin_validity', 'right_gaze_point_in_user_coordinate_system', 'right_gaze_point_validity', 'right_gaze_point_on_display_area'},...
	                		{'gaze_origin_3d', 'gaze_origin_validity', 'gaze_point_3d', 'gaze_point_validity', 'gaze_point_2d'});

	                end


	                % Remove invalid data from subsequent analyses
	                trialGazeData.gaze_origin_3d(trialGazeData.gaze_origin_validity==0,:) = NaN;
	                trialGazeData.gaze_point_3d(trialGazeData.gaze_point_validity==0,:) = NaN;
	                trialGazeData.gaze_point_2d(trialGazeData.gaze_point_validity==0,:) = NaN;

	                %savedEGdata(savedEGdata(:,8)==4, [2:7 9:10]) = NaN;
	                
	                % rescale 2d coordinates from proportion to pixels
	                trialGazeData.gaze_point_2d(:,1) = bsxfun(@times, trialGazeData.gaze_point_2d(:,1), res(1));
	                trialGazeData.gaze_point_2d(:,2) = bsxfun(@times, trialGazeData.gaze_point_2d(:,2), res(2));
	             	                
	                velSessionData(session).velPhaseData(phase).RawTrialData = trialGazeData.gaze_point_2d;
	                
	                % Interpolate Gaps in gaze position
	                GapsBeforeFill = sum(trialGazeData.gaze_point_validity==0)/size(trialGazeData.gaze_point_validity,1);
	                if GapsBeforeFill < 1 && GapsBeforeFill > 0 % valid data exists
	                    trialGazeData = fillMissing(trialGazeData, 75, sampleInterval); %interpolate gaps in data
	                end
	                GapsAfterFill = sum(trialGazeData.gaze_point_validity==0)/size(trialGazeData.gaze_point_validity,1);
	                
	                % apply moving average filter to smooth out noise in the gaze
	                % position samples
	                newX = zeros(size(trialGazeData.gaze_point_3d,1),1);
	                newY = zeros(size(trialGazeData.gaze_point_3d,1),1);
	                newZ = zeros(size(trialGazeData.gaze_point_3d,1),1);
	                
	                for ss = 1:size(trialGazeData.gaze_origin_3d,1)
	                    windowStep = 2; % average 2 samples either side of current data point
	                    windowCheck = 1;
	                    while windowCheck == 1
	                        windowStart = ss-windowStep;
	                        windowEnd = ss+windowStep;
	                        if windowStart < 1 || windowEnd > size(trialGazeData.gaze_point_3d,1) % if either the start or the end of the window aren't valid data points...
	                            windowStep = windowStep-1; % reduce the size of the moving average window
	                        elseif sum(trialGazeData.gaze_point_validity(windowStart:windowEnd)==0) > 1 % if any of the data within the moving window isn't valid (e.g. Blinks)
	                            windowStep = windowStep-1; % reduce the size of the moving average window
	                        else
	                            windowCheck = 0;
	                        end
	                    end
	                    
	                    if windowStep > -1 % if there is some valid data in the moving average window, average the X, Y, and Z gaze coordinates
	                        newX(ss) = mean(trialGazeData.gaze_point_3d(windowStart:windowEnd,1));
	                        newY(ss) = mean(trialGazeData.gaze_point_3d(windowStart:windowEnd,2));
	                        newZ(ss) = mean(trialGazeData.gaze_point_3d(windowStart:windowEnd,3));
	                    else
	                        newX(ss) = 0;
	                        newY(ss) = 0;
	                        newZ(ss) = 0;
	                    end
	                end
	                
	                % Update trialGazeData with the filtered gaze coords
	                trialGazeData.filtered_gaze_point = [newX, newY, newZ];

	                
	                % Start of I-VT classifier
	                velocity = zeros(1, size(trialGazeData.filtered_gaze_point,1)-sampleWindowLength); % set up velocity array
	                angle = zeros(1, size(trialGazeData.filtered_gaze_point,1)-sampleWindowLength); % set up angle array
	                movType = zeros(1, size(trialGazeData.filtered_gaze_point,1)-sampleWindowLength); % set up movement type array
	                
	                % check if gaze is within 100 px of scr centre for any sample in first 80ms (if discard set to on)
	                withinCentre = 1;
	                if discardOutsideFixation == 1
	                    withinCentre = 0;
	                    if size(trialGazeData.filtered_gaze_point,1) > 80/sampleInterval
	                        for ss = 1 : round(80/sampleInterval)
	                            if sqrt((trialGazeData.gaze_point_2d(ss,1)-fixationCoords(1))^2+(trialGazeData.gaze_point_2d(ss,2)-fixationCoords(2))^2) < 100 %pythagoras
	                                withinCentre = 1;
	                                break
	                            end
	                        end
	                    end
	                end
	                

	                aa = 1; %used as an index for stepping through the current eye movement
	                bb = 0; %used later on to keep track of different movements in trial
	                movList = {}; % this will hold the list of eye movements for this trial


	                for s = sampleWindowLength+1:size(trialGazeData.system_time_stamp,1)-sampleWindowLength % step through trial data with a moving window of ~20ms (8 samples)
	                    
	                    if sum(trialGazeData.gaze_point_validity(s-sampleWindowLength:s+sampleWindowLength)==4)>0 % if there is any missing data in the velocity window...
	                        velocity(s) = NaN;  %... velocity = NaN and move on
	                    else
	                        eyePoint = repmat(trialGazeData.gaze_origin_3d(s,:),2,1); % get XYZ eye point for sample. Doubled for vectorization later.
	                        gazePoint = trialGazeData.filtered_gaze_point([s-sampleWindowLength s+sampleWindowLength],:); % get filtered XYZ coords for start and end point of moving window
	                        gazeVector = eyePoint-gazePoint; % calculate vector from beginning and end gaze point to eye position in centre of time window
	                        angle(s) = acosd(dot(gazeVector(1,:),gazeVector(2,:))/(norm(gazeVector(1,:))*norm(gazeVector(2,:)))); %does fancy trig stuff to determine angle between vectors
	                        time = (trialGazeData.system_time_stamp(s+sampleWindowLength)-trialGazeData.system_time_stamp(s-sampleWindowLength))/1000000; %find time taken from beginning to end of time window in seconds
	                        velocity(s) = angle(s)/time; % determine angular velocity in degrees of visual angle/sec
	                    end
	                    
	                    	                    
	                    velSessionData(session).velPhaseData(phase).velTrialData{t} = velocity;
	                    velSessionData(session).velPhaseData(phase).TrialData{t} = trialGazeData.filtered_gaze_point;
	                    
	                    
	                    %I-VT classifier
	                    if velocity(s) > vThreshold
	                        movType(s) = 2; %mark as saccade
	                    elseif velocity(s) > -1
	                        movType(s) = 1; %mark as fixation
	                    else
	                        movType(s) = 0; %mark as gap
	                    end
	                    
	                    if movType(s) ~= movType(s-1) || s == size(trialGazeData.system_time_stamp,1)-sampleWindowLength % if detected movement is not the same as that on the previous sample, or we are at the end of the trial
	                        bb = bb + 1;
	                        if bb > 1 %add the previous eye movement to a higher level list
	                            endTime = (currentMov(end,1)+trialGazeData.system_time_stamp(s))/2; % figure out the end time of current movement in ms. Average of last sample in movement and next sample.
	                            movList(bb-1,:) = {movType(s-1) startTime endTime (endTime-startTime)/1000 latency currentMov currentMov(end,2:3)}; %[movement type, start time, end time, length of movement, movement details, x/y coords of last sample in movement]
	                            if movType(s-1) == 1   %if a fixation
	                                movList(bb-1,7) = {mean(currentMov(:,2:3),1)}; %calculate mean fixation point over entire fixation
	                            end
	                        end
	                        % collect info for new movement
	                        startTime = (trialGazeData.system_time_stamp(s)+trialGazeData.system_time_stamp(s-1))/2;
	                        latency = (startTime-trialGazeData.system_time_stamp(1))/1000;
	                        currentMov = [];
	                        currentMov(1,:) = [trialGazeData.system_time_stamp(s) trialGazeData.gaze_point_2d(s,:)];
	                        aa = 2;
	                    else %otherwise, add current sample info to current movement.
	                        currentMov(aa,:) = [trialGazeData.system_time_stamp(s) trialGazeData.gaze_point_2d(s,:)];
	                        aa = aa + 1;
	                    end
	                end


	                %% This is here so that you can view the filtered X and Y coordinates across time for each trial. Useful for debugging
	                    if graphVersion
	                        
	                        if ismember(t, trialsToCheck)
	                    
	                        msTime = (trialGazeData.system_time_stamp(:,1)-trialGazeData.system_time_stamp(1))/1000;
	                        acceleration = diff(velocity)/.0033;


	                        posFig = figure(1);
	                        subplot(2,2,1)
	                        plot(msTime, trialGazeData.gaze_point_2d(:,1), 'k')
	                        hold on
	                        plot(msTime, trialGazeData.gaze_point_2d(:,2), 'b')                   

	                        hold off
	                        axis([0, msTime(end), 0, 1920])
	                        title('X and Y Coordinates Across Time')
	                        xlabel('Time (ms)');
	                        ylabel('Pixel');
	                        legend('X Coordinates', 'Y Coordinates');

	                        subplot(2,2,3)
	                        plot(trialGazeData.gaze_point_2d(:,1), trialGazeData.gaze_point_2d(:,2), 'k')
	                        hold on
	                        for i = 1:length(stimCentre)
	                        	plot(stimCentre(i,1), stimCentre(i,2), 'kx')
	                        end

	                        targetTxt = '\leftarrow Target';
	                        distractTxt = '\leftarrow Distractor';
	                        text(stimCentre(targetLoc,1), stimCentre(targetLoc,2), targetTxt)
	                        text(stimCentre(distractLoc,1), stimCentre(distractLoc,2), distractTxt)
	                        hold off
	                        axis([0, 1920, 0, 1080])
	                        title('X and Y Coordinates')
	                        xlabel('X Coords')
	                        ylabel('Y Coords')

	                        subplot(2,2,2)
	                        plot(msTime(1:length(velocity)), velocity)
	                        hold on
	                        plot(msTime(1:length(velocity)), repmat(vThreshold, [length(velocity), 1]), '--k')
	                        hold off
	                        axis([0, msTime(end), 0, 100])
	                        title('Velocity Across Time')
	                        xlabel('Time (ms)')
	                        ylabel('Velocity (dva/s)')
	                        legend('Velocity', 'Velocity Threshold');
	                        
	                        subplot(2,2,4)
	                        plot(msTime(1:length(acceleration)), acceleration)
	                        hold on
	                        hold off
	                        axis([0, msTime(end), 0, 15000])
	                        title('Approximated Acceleration Across Time')
	                        xlabel('Time (ms)')
	                        ylabel('Acceleration (dva/s^2)')
	                        legend('Acceleration');
	                        pause;
	                        
	                        end
	                    
	                    end

	                if isempty(movList) == 0
	                    subMovements(subStep).sessionMovements(session).phaseMovements(phase).trialMovements(t) = {movList};
	                else
	                    subMovements(subStep).sessionMovements(session).phaseMovements(phase).trialMovements(t) = {0};
	                end
	                
	                
	                discardTrial = 0;
	                anticipatorySaccade = 0;
	                outsideFixation = 0;
	                noSaccades = 0;
	                noValidData = 0;
	                
	                %find first saccade
	                if isempty(movList) == 0
	                    saccList = cell2mat(movList(:,1)); % make movement list easier to work with
	                    saccIdx = find(saccList(:,1)==2); % index all saccades in movement list
	                    
	                    foundSaccade = 0;
	                    if isempty(saccIdx) == 0 % if there is at least one saccade detected
	                        aa = 1;
	                        idx = 0;
	                        while foundSaccade == 0 && aa <= sum(saccList(:,1)==2)
	                            idx = saccIdx(aa); % find index of the saccade
	                            foundSaccade = 1;
	                            saccadeLength = movList{idx,4}; % length of saccade in ms
	                            if saccadeLength < fixGapDurThresh; % if saccade is too short, likely a blink so disregard
	                                foundSaccade = 0;
	                                aa = aa + 1;
	                            end
	                            
	                        end

	                        saccadeLatency = movList{idx,5}; %find latency of first saccade
	                        
	                        if saccadeLatency < anticipationThreshold %if anticipatory saccade, or no samples near fixation within first 80ms, mark to be discarded
	                            discardTrial = 1;
	                            anticipatorySaccade = 1; %mark as anticipatory
	                        elseif withinCentre == 0
	                            discardTrial = 1;
	                            outsideFixation = 1;
	                        end
	                        
	                        
	                        endCoords = movList{idx,7}; %find end point of first saccade
	                        
	                        saccadeVector = [endCoords(1) - fixationCoords(1) endCoords(2) - fixationCoords(2)]; %calculate a vector from fixation point to saccade endpoint
	                        
	                        % determine the direction of the saccade
	                        direction = zeros(1,size(stimVector,1)+1);
	                        for aa = 1:size(stimVector,1);
	                            saccadeAngle = acosd(dot(saccadeVector,stimVector(aa,:))/(norm(saccadeVector)*norm(stimVector(aa,:)))); %calculate angle between fixation-stimulus vector and fixation-saccade endpoint vector
	                            if saccadeAngle <= withinAngle %if within a threshold angle
	                                direction(aa) = 1; %mark as going towards that stimulus
	                            end
	                        end
	                        
	                        if sum(direction) == 0 %if saccade is not within angle threshold of any ROI
	                            direction(end) = 1; % code as going "somewhere else"
	                        end
	                        
	                    end
	                    
	                    if foundSaccade == 0 %if no saccades detected in trial
	                        
	                        saccadeLatency = NaN;
	                        direction = ones(1,size(stimVector,1)+1)*99;
	                        discardTrial = 1;
	                        noSaccades = 1;
	                        
	                    end
	                    
	                else %if no valid eye data detected in trial
	                    saccadeLatency = NaN;
	                    direction = ones(1,size(stimVector,1)+1)*99;
	                    discardTrial = 1;
	                    noValidData = 1;
	                end
	                
	                saccadeSessionData(session).saccadePhaseData(phase).summarySaccadeData(t,:) = [sub, t, saccadeLatency, direction, discardTrial, anticipatorySaccade, outsideFixation, noSaccades, noValidData];
	                
	                movList = {};
	            end
	        end
	    
	    
	    
		    save([eyeDataFolder, '/processed_saccade_data/SummarySaccadeDataP',num2str(sub),'.mat'],'saccadeSessionData');
		    save([eyeDataFolder, '/raw_velocity_data/RawVelocityDataP', num2str(sub),'.mat'], 'velSessionData');
	    else
	    	disp(['Data for subject ', num2str(sub), ' session ', num2str(session), ' not found!'])
	    end
	end
end

save([eyeDataFolder, '/all_eye_movements/AllEyeMovements.mat'], 'subMovements');

end


function dataOut = fillMissing(dataIn, GapThresh, freqMS)
% dataIn = eyeGaze data to process (x,y,validity)
% GapThresh = time window in ms for acceptable gaps
% freqMS = ms value of each timestamp


% Interpolation of gaps
while dataIn.gaze_point_validity(1) == 0 % remove gaps at start
    dataIn.gaze_point_validity(1) = [];
    dataIn.system_time_stamp(1) = [];
    dataIn.gaze_point_2d(1,:) = [];
    dataIn.gaze_point_3d(1,:) = [];
    dataIn.gaze_origin_3d(1,:) = [];
end
while dataIn.gaze_point_validity(1) == 0 % remove gaps at end
    dataIn.gaze_point_validity(end) = [];
    dataIn.system_time_stamp(end) = [];
    dataIn.gaze_point_2d(end,:) = [];
    dataIn.gaze_point_3d(end,:) = [];
    dataIn.gaze_origin_3d(end,:) = [];
end

% this section works out what positions needs to be filled and provides
% start and end points for each fill, stored in iFills
iFills = zeros(size(dataIn.gaze_point_validity,1),2);
intCnt = 0;
checkPos = 2;
while checkPos < size(dataIn.gaze_point_validity,1) % check each position in turn
    endPos = checkPos; % set end to current check
    if dataIn.gaze_point_validity(checkPos)==0 % if missing (otherwise increase check position)
        while dataIn.gaze_point_validity(endPos)==0 % step through until valid data is found
            endPos = endPos + 1; % increase end position
        end
        if endPos-checkPos < (GapThresh/freqMS) % if that gap is smalle enough
            intCnt = intCnt + 1; % this is a new fill
            iFills(intCnt,:) = [checkPos-1 endPos]; % add details of fill to the array
        end
        checkPos = endPos + 1; % go to next check position beyond the end
    else
        checkPos = checkPos + 1;
    end
end
iFills(intCnt+1:end,:) = []; % remove empty rows of array

% use values to interpolate
for r = 1:size(iFills,1)
    gaze_origin_3d = [];
    gaze_point_3d = [];
    gaze_point_2d = [];
    intSteps = 0:1/(iFills(r,2)-iFills(r,1)):1; % calculate appropriate distribution across fill range
    
    gaze_origin_3d(:,1) = dataIn.gaze_origin_3d(iFills(r,1),1) + (dataIn.gaze_origin_3d(iFills(r,2),1)-dataIn.gaze_origin_3d(iFills(r,1),1))*intSteps; % interpolation of x
    gaze_origin_3d(:,2) = dataIn.gaze_origin_3d(iFills(r,1),2) + (dataIn.gaze_origin_3d(iFills(r,2),2)-dataIn.gaze_origin_3d(iFills(r,1),2))*intSteps; % interpolation of y
    gaze_origin_3d(:,3) = dataIn.gaze_origin_3d(iFills(r,1),3) + (dataIn.gaze_origin_3d(iFills(r,2),3)-dataIn.gaze_origin_3d(iFills(r,1),3))*intSteps; % interpolation of z\
    gaze_point_3d(:,1) = dataIn.gaze_point_3d(iFills(r,1),1) + (dataIn.gaze_point_3d(iFills(r,2),1)-dataIn.gaze_point_3d(iFills(r,1),1))*intSteps;
    gaze_point_3d(:,2) = dataIn.gaze_point_3d(iFills(r,1),2) + (dataIn.gaze_point_3d(iFills(r,2),2)-dataIn.gaze_point_3d(iFills(r,1),2))*intSteps;
    gaze_point_3d(:,3) = dataIn.gaze_point_3d(iFills(r,1),3) + (dataIn.gaze_point_3d(iFills(r,2),3)-dataIn.gaze_point_3d(iFills(r,1),3))*intSteps;
    gaze_point_2d(:,1) = dataIn.gaze_point_2d(iFills(r,1),1) + (dataIn.gaze_point_2d(iFills(r,2),1)-dataIn.gaze_point_2d(iFills(r,1),1))*intSteps;
    gaze_point_2d(:,2) = dataIn.gaze_point_2d(iFills(r,1),2) + (dataIn.gaze_point_2d(iFills(r,2),2)-dataIn.gaze_point_2d(iFills(r,1),2))*intSteps;
    dataIn.gaze_origin_3d(iFills(r,1):iFills(r,2),:) = [gaze_origin_3d];
    dataIn.gaze_point_3d(iFills(r,1):iFills(r,2),:) = [gaze_point_3d];
    dataIn.gaze_point_2d(iFills(r,1):iFills(r,2),:) = [gaze_point_2d];
end

dataOut = dataIn;

end
