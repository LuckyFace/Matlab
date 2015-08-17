function Event2Mat(sessionFolder)
% Event2Mat Converts data from Neuralynx NEV files to Matlab mat files

% 주어진 sessionFolder 안에 있는 nev file의 목록을 만든다
narginchk(0, 1);
if nargin == 0
    eventFiles = FindFiles('*.nev','CheckSubdirs',0);
elseif nargin == 1
    if ~iscell(sessionFolder)
        disp('Input argument is wrong. It should be cell array.');
        return;
    elseif isempty(sessionFolder)
        eventFiles = FindFiles('*.nev','CheckSubdirs',0);
    else
        nFolder = length(sessionFolder);
        eventFiles = cell(0,1);
        for iFolder = 1:nFolder
            if exist(sessionFolder{iFolder},'dir')
                cd(sessionFolder{iFolder});
                eventFiles = [eventFiles;FindFiles('Events.nev','CheckSubdirs',1)];
            end
        end
    end
end
if isempty(eventFiles)
    disp('Event file does not exist!');
    return;
end

nFile = length(eventFiles);
for iFile = 1:nFile
    disp(['Analyzing ',eventFiles{iFile}]);
    cd(fileparts(eventFiles{iFile}));
    
    [timeStamp, eventString] = Nlx2MatEV(eventFiles{iFile}, [1 0 0 0 1], 0, 1, []);
    timeStamp = timeStamp';
    
    % Task
    trialOnsetIndex = strncmp(eventString, 'Cue', 3);
    trialOffsetIndex = strcmp(eventString, 'Baseline');
    rewardIndex = strcmp(eventString, 'Reward') | strcmp(eventString, 'Non-reward');
    offsetIndex = reshape(find(strncmp(eventString, 'TTL Input on AcqSystem1_0 board 0 port 2', 40)), 3, []);
    lickOnsetIndex = strcmp(eventString, 'Sensor');
    lickOffsetIndex = strncmp(eventString, 'TTL Input on AcqSystem1_0 board 0 port 3', 40);
    
    nTrial = sum(trialOnsetIndex);
    
    trialOnsetTime = timeStamp(trialOnsetIndex);
    cueOnsetTime = timeStamp(offsetIndex(1,:));
    delayOnsetTime = timeStamp(offsetIndex(2,:));
    rewardOnsetTime = timeStamp(rewardIndex);
    rewardOffsetTime = timeStamp(offsetIndex(3,:));
    lickOnsetTime = timeStamp(lickOnsetIndex);
    lickOffsetTime = timeStamp(lickOffsetIndex);
    trialTime = timeStamp(trialOffsetIndex);
    
    cue = cellfun(@(x) str2num(x(4)), eventString(trialOnsetIndex), 'UniformOutput', true);
    reward = strcmp(eventString(rewardIndex), 'Reward');
    
    
    cueA = cue==1;
    cueB = cue==2;
    cueC = cue==3;
    cueD = cue==4;
        
    cueA_reward = (cue==1) & (reward==1);
    cueB_reward = (cue==2) & (reward==1);
    cueC_reward = (cue==3) & (reward==1);
    cueD_reward = (cue==4) & (reward==1);
    cueA_nonreward = (cue==1) & (reward==0);
    cueB_nonreward = (cue==2) & (reward==0);
    cueC_nonreward = (cue==3) & (reward==0);
    cueD_nonreward = (cue==4) & (reward==0);
    
    cueIndex = [cueA, cueB, cueC, cueD];
    cueResult = sum(cueIndex);
    trialIndex = [cueA_reward, cueA_nonreward, cueB_reward, cueB_nonreward, ...
        cueC_reward, cueC_nonreward, cueD_reward, cueD_nonreward];
    eventTime = [trialOnsetTime, cueOnsetTime, delayOnsetTime, rewardOnsetTime];
    trialResult = sum(trialIndex);
    
    binSize = 0.01; % unit: second;
    lickBin = (0:binSize:8);
    lickHist = zeros(nTrial,length(lickBin));
    resolution = 10;
    for iTrial = 1:nTrial
        lickHist(iTrial,:) = histc(lickOnsetTime, trialTime(iTrial) + lickBin*1000000)/binSize;
    end
    
    lickMeanConv = zeros(4, length(lickBin));
    lickPSTHSem = zeros(4, length(lickBin));
    for iType = 1:4
        lickMean = mean(lickHist(cueIndex(:,iType),:));
        lickSem = std(lickHist(cueIndex(:,iType),:))/sqrt(cueResult(iType));
        lickMeanConv(iType,:) = conv(lickMean, fspecial('Gaussian', [1 5*resolution],resolution), 'same');
        lickSemConv(iType,:) = conv(lickSem, fspecial('Gaussian', [1 5*resolution],resolution), 'same');
    end
    
    lickSemBin = [lickBin flip(lickBin)];
    lickSemConv = [lickMeanConv-lickSemConv flip(lickMeanConv+lickSemConv,2)];
    
    save('Events.mat', ...
        'nTrial', 'trialIndex', 'eventTime', 'trialTime', 'trialResult', ...
        'cueIndex', 'cueResult', ...
        'cue', 'reward', 'lickOnsetTime','lickOffsetTime', ...
        'lickBin','lickMeanConv','lickSemBin','lickSemConv');   
end
disp('Done!');