function psthCC(sessionFolder)
% psthCC Converts data from MClust t files to Matlab mat files

% Task variables
win = [-1 9]*10^3; % unit: msec, window for binning
window = [0 8]; % unit: sec, final view window
binSize = 0.01; % unit: sec, = 10 msec
resolution = 10; % sigma = resolution * binSize = 100 msec

% Tag variables
winTag = [-20 100]; % unit: msec
binTagSize = 2;
binTag = winTag(1):binTagSize:winTag(2);

% Find files
switch nargin
    case 0 % Input이 없는 경우 그냥 폴더안의 t 파일을 검색
        ttFile = FindFiles('T*.t','CheckSubdirs',0); % Subfolder는 검색하지 않는다
    case 1 % Input이 있는 경우
        if ~iscell(sessionFolder) % 셀 array인지 확인
            disp('Input argument is wrong. It should be cell array.');
            return;
        elseif isempty(sessionFolder) % Cell이 맞지만 텅 비었으면 그냥 폴더 안의 t파일 검색
            ttFile = FindFiles('T*.t','CheckSubdirs',1);
        else % Cell이 맞고 빈 array가 아니면, 차례대로 cell 내용물 확인
            nFolder = length(sessionFolder);
            ttFile = cell(0,1);
            for iFolder = 1:nFolder
                if exist(sessionFolder{iFolder})==7 % 폴더이면 그 아래 폴더들의 t파일 검색
                    cd(sessionFolder{iFolder});
                    ttFile = [ttFile;FindFiles('T*.t','CheckSubdirs',1)];
                elseif strcmp(sessionFolder{iFolder}(end-1:end),'.t') % t파일이면 바로 합친다.
                    ttFile = [ttFile;sessionFolder{iFolder}];
                end
            end
        end
end
if isempty(ttFile)
    disp('TT file does not exist!');
    return;
end
ttData = LoadSpikes(ttFile,'tsflag','ts','verbose',0);

nCell = length(ttFile);
for iCell = 1:nCell
    disp(['### Analyzing ',ttFile{iCell},'...']);
    [cellPath,cellName,~] = fileparts(ttFile{iCell});
    cd(cellPath);

    % Event variables
    load('Events.mat');
    yTemp = [0:nTrial-1; 1:nTrial; NaN(1,nTrial)]; % ypt 만들 때 사용
    resultSum = [0 cumsum(trialResult)];
    
    spikeData = Data(ttData{iCell})/10; % unit: msec
    spikeTime = cell(nTrial,1);
    
    % Firing rate
    fr_base = sum(histc(spikeData,baseTime))/diff(baseTime/1000);
    fr_task = sum(histc(spikeData,taskTime))/diff(taskTime/1000);
        
    % Making raster, PSTH
    for iTrial = 1:nTrial
        timeIndex = [];
        if isnan(eventTime(iTrial,1)); continue; end;
        [~,timeIndex] = histc(spikeData,eventTime(iTrial,1)+win);
        if isempty(timeIndex); continue; end;
        spikeTime{iTrial} = spikeData(logical(timeIndex))-eventTime(iTrial,1);
    end

    % Making raster points
    spikeBin = win(1)/10^3:binSize:win(2)/10^3; % unit: sec
    totalHist = histc(cell2mat(spikeTime),spikeBin)/(binSize*nTrial);
    fireMean = mean(totalHist);
    fireStd = std(totalHist);
    
    nCue = length(trialResult);
    nSpikeBin = length(spikeBin);
    xpt = cell(1,nCue);
    ypt = cell(1,nCue);
    spikeHist = zeros(nCue,nSpikeBin);
    spikeConv = zeros(nCue,nSpikeBin);

    for iCue = 1:nCue % Cue A_rewarded (Ay), An, By, Bn, Cy, Cn, Dy, Dn
        nSpikePerTrial = cellfun(@length,spikeTime(trialIndex(:,iCue)));
        nSpikeTotal = sum(nSpikePerTrial);
        if trialResult(iCue)==0 || nSpikeTotal == 0
            xpt{iCue} = [];
            ypt{iCue} = [];
            continue;
        end

        spikeTemp = cell2mat(spikeTime(trialIndex(:,iCue)))';
        xptTemp = [spikeTemp;spikeTemp;NaN(1,nSpikeTotal)];
        xptTemp = xptTemp(:);

        yptTemp = [];
        for iy = 1:trialResult(iCue)
            yptTemp = [yptTemp repmat(yTemp(:,resultSum(iCue)+iy),1,nSpikePerTrial(iy))];
        end
        yptTemp = yptTemp(:);
        xpt{iCue} = xptTemp/1000;
        ypt{iCue} = yptTemp;

        % Making PSTH
        spkhist_temp = histc(spikeTemp/1000,spikeBin)/(binSize*trialResult(iCue));
        spkconv_temp = conv(spkhist_temp,fspecial('Gaussian',[1 5*resolution],resolution),'same');
        spikeHist(iCue,:) = spkhist_temp;
        spikeConv(iCue,:) = spkconv_temp;
    end
    peth = spikeHist;
    pethconv = spikeConv;
    zpeth = (spikeHist-fireMean)/fireStd;
    zpethconv = (spikeConv-fireMean)/fireStd;

    save([cellName,'.mat'],...
        'fr_base','fr_task',...
        'xpt','ypt',...
        'peth','pethconv',...
        'zpeth','zpethconv',...
        'window','spikeBin');
    
    %% Tagging
    nPulse = length(blueOnsetTime);
    yTagTemp =[0:(nPulse-1);1:nPulse;NaN(1,nPulse)];
    spkTemp = cell(nPulse,1);
    
    for iPulse = 1:nPulse
        timeIndex = [];
        [~,timeIndex] = histc(spikeData,blueOnsetTime(iPulse)+winTag);
        if isempty(timeIndex); continue; end;
        spkTemp{iPulse} = spikeData(logical(timeIndex)) - blueOnsetTime(iPulse);
    end
    spikeTagTime = cell2mat(spkTemp)';
    tagHist = histc(spikeTagTime,binTag)/(binTagSize/1000*nPulse);
    nTagPerTrial = cellfun(@length,spkTemp);
    nTagTotal = sum(nTagPerTrial);
    xpttag = []; ypttag = [];
    if nTagTotal~=0;
        xpttag = [spikeTagTime;spikeTagTime;NaN(1,nTagTotal)];
        xpttag = xpttag(:);
        for iytag = 1:nPulse
            ypttag = [ypttag repmat(yTagTemp(:,iytag),1,nTagPerTrial(iytag))];
        end
        ypttag = ypttag(:);
    else
        xpttag = []; ypttag = [];
        tagHist = zeros(1,length(binTag));      
    end
    save([cellName,'.mat'],...
        'xpttag','ypttag','binTag','tagHist','-append');
end
disp('### Making Raster, PSTH is done!');