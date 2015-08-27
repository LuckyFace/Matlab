function psthplotCC(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Object: Mat으로 저장된 정보로 PSTH/Raster/Tagging 등의 그림을 그린다.
% Author: Dohoung Kim
% First written: 2014/09/19
%
% Ver 2.0 (2015. 1. 7)
%   Function 수정
%   Figure file은 처음 실행한 directory에서 저장되도록
%   Modulation하지 않은 파일만 분석한다. 
% Ver 2.1 (2015. 1. 9)
%   Log-rank test plot 추가
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables
win = [-20 100]; % Tagging period
p_value = 0.05; % Threshold for tagging

% Plot properties
lineclr={[0.8 0 0], [0.8 0.4 0], ...
        [0 0 0.8], [0 0.4 0.8], ...
        [1 0.8125 0], [1 0.8125 0], ...
        [0 0.8125 1], [0 0.8125 1]};
linewth=[2 1 2 1 2 1 2 1];
linestl={'-','-','-','-','-','-','-','-'};

% Find files
switch nargin
    case 0 % Input이 없는 경우 그냥 폴더안의 mat 파일을 검색
        matFile = FindFiles('T*.mat','CheckSubdirs',0); % Subfolder는 검색하지 않는다
    case 1 % Input이 있는 경우
        if ~iscell(folders) % 셀 array인지 확인
            disp('Input argument is wrong. It should be cell array.');
            return;
        elseif isempty(folders) % Cell이 맞지만 텅 비었으면 그냥 폴더 안의 mat파일 검색
            matFile = FindFiles('T*.mat','CheckSubdirs',1); % Subfolder 검색.
        else % Cell이 맞고 빈 array가 아니면, 차례대로 cell 내용물 확인
            nfolder = length(folders);
            matFile = cell(0,1);
            for ifolder = 1:nfolder
                if exist(folders{ifolder})==7 % 폴더이면 그 아래 폴더들의 mat파일 검색
                    cd(folders{ifolder});
                    matFile = [matFile;FindFiles('T*.mat','CheckSubdirs',1)];
                elseif strcmp(folders{ifolder}(end-3:end),'.mat') % mat파일이면 바로 합친다.
                    matFile = [matFile;folders{ifolder}];
                end
            end
        end
end
if isempty(matFile)
    disp('Mat file does not exist!');
    return;
end
nFiles = length(matFile);
rtdir = pwd;

% Plot position
% fHandle = figure('PaperUnits','centimeters','PaperPosition',[2 2 8.9 6.88]);
nRow = 5;
startPoint = [0.5 0.1];
figWidth = [0.4 0.8];
interval_y = 0.05;
interFig_y = 0.01;
dx = figWidth(1);
dy = (figWidth(2)-(nRow-1)*interval_y-nRow*interFig_y)/(2*nRow);
axpt = zeros(nRow,2,4);
for iRow = 1:nRow
    axpt(iRow,1,:) = [startPoint(1),...
        startPoint(2)+(nRow-iRow)*(2*dy+interFig_y+interval_y)+interFig_y+interval_y,...
        dx dy];
    axpt(iRow,2,:) = [startPoint(1),...
        startPoint(2)+(nRow-iRow)*(2*dy+interFig_y+interval_y),...
        dx dy];
end

load(matFile{1},'spikeBin','window');
iPage = 0;
for iFiles = 1:nFiles
    [cellcd,cellname,~] = fileparts(matFile{iFiles});
    
    cd(cellcd);
    
    load('Events.mat','nTrial','trialResult','blueOnsetTime');
    lightTrial = length(blueOnsetTime);
    load(matFile{iFiles},'xpt','ypt','pethconv');
    ylims = 0;
    
%% PETH and Raster for epoch
% Raster
    axes('Position',axpt(mod(iFiles-1,nRow)+1,1,:)); % Raster
    hold on;
    plot([0.5 0.5],[0.01 nTrial],'LineStyle','-','LineWidth',0.3,'Color',[0.8 0.8 0.8]);
    plot([1.5 1.5],[0.01 nTrial],'LineStyle','-','LineWidth',0.3,'Color',[0.8 0.8 0.8]);
    plot([4 4],[0.01 nTrial],'LineStyle','-','LineWidth',0.3,'Color',[0.8 0.8 0.8]);
    for iChoice = find(trialResult~=0)
        plot(xpt{iChoice},ypt{iChoice},...
            'LineStyle','-','LineWidth',0.2,'Color',lineclr{iChoice});
    end
    set(gca,'box','off','TickDir','out','LineWidth',0.2,'FontSize',4, ...
        'XLim',window,'XTick',[], ...
        'YLim',[0 nTrial],'YTick',[0 nTrial],'YTickLabel',{'      ',nTrial});
    ylabel('Trial','FontSize',4);

% PETH
    axes('Position',axpt(mod(iFiles-1,nRow)+1,2,:)); % PSTH
    hold on;
    plot([0.5 0.5],[0.01 nTrial],'LineStyle','-','LineWidth',0.3,'Color',[0.8 0.8 0.8]);
    plot([1.5 1.5],[0.01 nTrial],'LineStyle','-','LineWidth',0.3,'Color',[0.8 0.8 0.8]);
    plot([4 4],[0.01 nTrial],'LineStyle','-','LineWidth',0.3,'Color',[0.8 0.8 0.8]);
    for jChoice = find(trialResult~=0)
        plot(spikeBin,pethconv(jChoice,:),...
            'LineStyle',linestl{jChoice},'LineWidth',linewth(jChoice),'Color',lineclr{jChoice});
        
    end
    ylims = ceil(max(pethconv(:))*1.1);
    set(gca,'box','off','TickDir','out','LineWidth',0.2,'FontSize',4, ...
        'XLim',window, 'XTick',[0 0.5 1.5 4 8], ...
        'YLim',[0 ylims], 'YTick',[0 ylims], 'YTickLabel',{'     ',ylims});
    if (mod(iFiles,nRow)==0 | iFiles == nFiles)
                xlabel('Time (s)','FontSize',5);
            end
    ylabel('Rate (Hz)','FontSize',4);

% File name
axes('Position',[axpt(1,2,1) axpt(mod(iFiles-1,nRow)+1,2,2)-0.03 0.4 0.1])
text(0,0,matFile{iFiles},'FontSize',4,'interpreter','none');
set(gca,'Visible','off');

%% Tagging
    clear xpttag ypttag bintag taghist time_tagstat H1_tagstat H2_tagstat p_tagstat;
    load(matFile{iFiles},'xpttag','ypttag','binTag','tagHist','time_tagstat','H1_tagstat','H2_tagstat','p_tagstat');
    if ~isempty(xpttag)
        axes('Position',[axpt(mod(iFiles-1,nRow)+1,1,1)-0.35 axpt(mod(iFiles-1,nRow)+1,1,2) dx/5 dy]);
            hold on;
            plot(xpttag,ypttag,...
                    'Marker','.','MarkerSize',2,'Color','k');
            set(gca,'box','off','TickDir','out','LineWidth',0.2,'FontSize',4, ...
                'XLim',win,'XTick',[], ...
                'YLim',[0 lightTrial],'YTick',[0 lightTrial],'YTickLabel',{[],lightTrial});
            ylabel('Trials');

        axes('Position',[axpt(mod(iFiles-1,nRow)+1,2,1)-0.35 axpt(mod(iFiles-1,nRow)+1,2,2) dx/5 dy]);
            hold on;
            h = bar(binTag,tagHist,'histc');
            ylims2 = ceil(max(tagHist(:))*1.1);
            if ylims2 ==0; ylims2=1; end;
            set(h,'FaceColor','k','EdgeAlpha',0);
            set(gca,'XLim',win);
            set(gca,'YLim',[0 ylims2],'YTick',[0 ylims2],'YColor','k');
            set(gca,'Box','off','TickDir','out','LineWidth',0.2,'FontSize',4, ...
                'XLim',win, 'XTick',[-20 0 100]);
            if (mod(iFiles,nRow)==0 | iFiles == nFiles)
                xlabel('Time (ms)');
            end
            ylabel('Rate (Hz)','FontSize',4);

        axes('Position',[axpt(mod(iFiles-1,nRow)+1,1,1)-0.22 axpt(mod(iFiles-1,nRow)+1,2,2) dx/5 dy]);
            hold on;
            stairs(time_tagstat,H1_tagstat,'LineStyle','-','LineWidth',0.5,'Color',[0 0.66 1]);
            stairs(time_tagstat,H2_tagstat,'LineStyle',':','LineWidth',0.5,'Color',[0 0 0]);
            ylimst = max([H1_tagstat;H2_tagstat])*1.1;
            if isempty(ylimst) | ylimst==0; ylimst=1; end;
            set(gca,'Box','off','TickDir','out','FontSize',4);
            set(gca,'XLim',[0 10],'XTick',[0 10]);
            set(gca,'YLim',[0 ylimst],'YTick',[0 ylimst],'YTickLabel',{[0],num2str(ylimst,1)});
            if (mod(iFiles,nRow)==0 | iFiles == nFiles)
                xlabel('Time (ms)');
            end
            ylabel('H(t)');

        axes('Position',[axpt(mod(iFiles-1,nRow)+1,1,1)-0.22 axpt(mod(iFiles-1,nRow)+1,1,2) dx/5 dy]);
            hold on;
            text(0,0.5,['p value for tag test: ',num2str(p_tagstat,2)],'FontSize',4,'interpreter','none');
            if p_tagstat <= p_value
                if H1_tagstat(end) >= H2_tagstat(end)
                    text(0,0.3,['Activated neuron'],'FontSize',4,'interpreter','none');
                else
                    text(0,0.3,['Inhibited neuron'],'FontSize',4,'interpreter','none');
                end
            end
            set(gca,'visible','off');
    end
    
%% Waveform
     load(matFile{iFiles},'spkwv','hfvwth','spkpvr','fr_base','fr_task');
     ylims3 = [min(spkwv(:)) max(spkwv(:))];
    
    for ich = 1:4
        axes('Position',[axpt(mod(iFiles-1,nRow)+1,1,1)-0.13+(ich-1)*dx/24*0.8 axpt(mod(iFiles-1,nRow)+1,2,2) dx/24*0.8 dy*0.8]);
        plot(spkwv(ich,:),'LineStyle','-','LineWidth',0.1,'Color',[0.2 0.2 0.2]);
        
        set(gca,'XLim',[1 32]);
        set(gca,'YLim',ylims3);
        set(gca,'visible','off');
        if ich==4
            line([24 32],[ylims3(2)-50 ylims3(2)-50],'color','k','LineWidth',0.2); 
            line([24 24],[ylims3(2)-50 ylims3(2)],'color','k','LineWidth',0.2);
        end
    end
    
     axes('Position',[axpt(mod(iFiles-1,nRow)+1,1,1)-0.13 axpt(mod(iFiles-1,nRow)+1,1,2) dx/6 dy]);
        hold on;
         text(0,0.9,['Base firing rate: ',num2str(fr_base,3)],'FontSize',4,'interpreter','none');
         text(0,0.7,['Task firing rate: ',num2str(fr_task,3)],'FontSize',4,'interpreter','none');
         text(0,0.5,['Half-valley width: ',num2str(hfvwth,3),' us'],'FontSize',4,'interpreter','none');
         text(0,0.3,['Peak valley ratio: ',num2str(spkpvr,3)],'FontSize',4,'interpreter','none');
         set(gca,'visible','off');
       
%% Title and save image        
    if (mod(iFiles,nRow)==0 | iFiles == nFiles)
        axes('Position',[axpt(1,1,1) axpt(1,1,2)+dy+0.006 dx 0.05]);
            text(0.25,0,'Preparation','FontSize',4,'interpreter','none','HorizontalAlignment','center');
            text(1,0,'Cue','FontSize',4,'interpreter','none','HorizontalAlignment','center');
            text(2.75,0,'Delay','FontSize',4,'interpreter','none','HorizontalAlignment','center');
            text(4.5,0,'Reward','FontSize',4,'interpreter','none','HorizontalAlignment','center');
            set(gca,'Visible','off','XLim',[0 8]);
            
        cd(rtdir);
        if nargin == 0
            cd('..');
        end
        iPage = iPage + 1;
        cell_filename = regexp(cellcd,'\','split');
        cellfile = strcat(cell_filename(end-1),'_',cell_filename(end),'_Raster_',num2str(iPage),'.tif');
        print(gcf,'-dtiff','-r600',cellfile{1});
        clf;
    end
end
close all;
cd(rtdir);