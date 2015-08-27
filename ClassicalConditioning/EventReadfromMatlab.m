clc; clear all; close all;

cd('D:\Cloud\project\classical_conditioning\data\PVMD1\');

matFile = FindFiles('PV*.mat');
nFile = length(matFile);

lineColor = {[1 0 0], [0 0 1]};
rewardStyle = {'-','--'};

for iFile = nFile
    load(matFile{iFile});
    
    if stateTime(1,1) == 0 || isnan(stateTime(1,1))
        stateTime(1,1) = stateTime(1,2) - 500006;
    end
    
    binSize = 0.01;
    bin = 0:binSize:8;
    nBin = length(bin);
    startTime = (repmat(stateTime(:,1),1,nBin) + repmat(bin*1000000,nTrial,1))';
    
    binCount = histc(lickTime(:,1),startTime);
    binCount = reshape(binCount,nBin,nTrial)';
    
    PSTH = zeros(4,nBin);
    PSTHconv = zeros(4,nBin);
    resolution = 10;
    hold on;
    for iCue = 1:2
        for iReward = 1:2
            index = 2*(iCue-1)+iReward;
            PSTH(index,:) = mean(binCount(odorCue==(iCue-1) & waterReward==(2-iReward),:))/binSize;
            PSTHconv(index,:) = conv(PSTH(index,:), fspecial('Gaussian', [1 5*resolution], resolution), 'same');
            
            plot(bin,PSTHconv(index,:),'Color',lineColor{iCue},'LineStyle',rewardStyle{iReward});
        end
    end
    
end
        