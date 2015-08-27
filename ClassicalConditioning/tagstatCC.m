function tagstat(folders)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Object: Tagging 시 통계 검정
% Author: Dohoung Kim
% First written: 2014/09/19
%
% Ver 2.0 (2015. 1. 1)
%   Tagging psth / 검정도 같이 통합
% Ver 2.1 (2015. 1. 8)
%   Log-rank test로 변경
% Ver 2.2 (2015. 1. 15)
%   200 trial 마다 선택해서 분석한 결과도 넣기
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables
binrange = 10;
binnum = 40;
slide = 1;

bin = [-binnum*binrange:binrange:0];
bin(end) = slide;
nbin = length(bin);

% Find files
switch nargin
    case 0 % Input이 없는 경우 그냥 폴더안의 t 파일을 검색
        ttfile = FindFiles('T*.t','CheckSubdirs',0); % Subfolder는 검색하지 않는다
        if isempty(ttfile)
            disp('TT file does not exist!');
            return;
        end
    case 1 % Input이 있는 경우
        if ~iscell(folders) % 셀 array인지 확인
            disp('Input argument is wrong. It should be cell array.');
            return;
        elseif isempty(folders) % Cell이 맞지만 텅 비었으면 그냥 폴더 안의 t파일 검색
            ttfile = FindFiles('T*.t','CheckSubdirs',1);
            if isempty(ttfile)
                disp('Cluster does not exist!');
                return;
            end
        else % Cell이 맞고 빈 array가 아니면, 차례대로 cell 내용물 확인
            nfolder = length(folders);
            ttfile = cell(0,1);
            for ifolder = 1:nfolder
                if exist(folders{ifolder})==7 % 폴더이면 그 아래 폴더들의 t파일 검색
                    cd(folders{ifolder});
                    ttfile = [ttfile;FindFiles('T*.t','CheckSubdirs',1)];
                elseif strcmp(folders{ifolder}(end-1:end),'.t') % t파일이면 바로 합친다.
                    ttfile = [ttfile;folders{ifolder}];
                end
            end
            if isempty(ttfile)
                disp('TT file does not exist!');
                return;
            end
        end
end
ncell = length(ttfile);
ttdata = LoadSpikes(ttfile,'tsflag','ts');

for icell = 1:ncell
    [cellpath,cellname,~] = fileparts(ttfile{icell});
    cd(cellpath);
    disp(['### Log-rank test: ',ttfile{icell}]);
    
    load('Events.mat','blueOnsetTime');
    
    outbin = find(diff(blueOnsetTime)<=400);
    outbin = [outbin;(outbin+1)];
    blueOnsetTime(outbin(:))=[];
    nlight = length(blueOnsetTime);
    
    binmat = ones(nlight,nbin)*diag(bin);
    lightbin = (repmat(blueOnsetTime',nbin,1)+binmat');
    time = zeros(nbin,nlight);
    censor = zeros(nbin,nlight);
    
    spkdata = Data(ttdata{icell})/10;

    for ilight=1:nlight
        for ibin=1:nbin
            idx=find(spkdata>lightbin(ibin,ilight),1,'first');
            if isempty(idx)
                time(ibin,ilight)=binrange;
                censor(ibin,ilight)=1;
            else
                time(ibin,ilight)=spkdata(idx)-lightbin(ibin,ilight);
                if time(ibin,ilight)>binrange
                    time(ibin,ilight)=binrange;
                    censor(ibin,ilight)=1;
                end
            end     
        end
    end
    base = [reshape(time(1:binnum,:),1,[]);reshape(censor(1:binnum,:),1,[])]';
    test = [time(end,:);censor(end,:)]';
    
    [p_tagstat,time_tagstat,H1_tagstat,H2_tagstat] = logrank(test,base);

    save([cellname,'.mat'],...
        'p_tagstat','time_tagstat','H1_tagstat','H2_tagstat',...
        '-append');
end
disp('### Log-rank test done!');