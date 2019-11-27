%% Download Data 

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project
%first, check that work environment is set, if not, creat the folders
DataOutFolder ='./HDW/2010' ;
if isempty(ls(DataOutFolder))
    mkdir(DataOutFolder);
%     mkdir([DataOutFolder '/Unzipped']);
end

StartYear = 2006; %start year and month of dataset
% StartMonth = 9; 

EndYear = 2015;%end year and month of dataset
% EndMonth= 12;

% setting of start/end year/month can be improved by obtaining from input
str = strings(10,366);
i=0;
for Year = StartYear:EndYear %loop over years
    i = i+1;
    j=0;
    Month1 = 1;
    Month12 = 12;
    
    for Month = Month1:Month12
        % number of days in the month
        if(Month==1)
            Day1 = 4;
        else
            Day1 =1;
        end
        
        if (Month ==12)
            Day2 = 29;
        elseif (Month==2)
            if mod(Year,4) == 0 %in the first year, start in the start month rather than January
                Day2 = 29;
            else
                Day2=28;
            end
        elseif ((mod(Month,2)==0 && Month>7 && Month<12) ||(mod(Month,2)==1 && Month<=7))
            Day2 = 31;
        else
            Day2 = 30;
        end
        
        % download daily file
        for Day = Day1:Day2
            % for cfwi
%           FileName = ['cfwi_' num2str(Year) '_' num2str(Month,'%02d') ...
%                '_' num2str(Day,'%02d') '.nc']; %construct the file name
            % for hdw
            FileName = ['out_' num2str(Year) num2str(Month,'%02d') num2str(Day,'%02d') '.nc'];
            
            j = j+1;
            disp(FileName);
            str(i,j)= FileName;
            % for cfwi
%           outfilename = websave([DataOutFolder '/' FileName],['http://eams3.usfs.msu.edu/study/2019-WRF-Land-Use-pyth/cfwi2/future/' num2str(Year) '/' FileName]);
            % for HDW
            outfilename = websave([DataOutFolder '/' FileName],['http://eams3.usfs.msu.edu/study/2019-HDW_WRFLULC/LULCC-2011/' FileName]);
        end        
    end  
end

%% Group files by year & month
% HDW-Current
FolderPath = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\HDW\2010';
Index = 'HDW';
ReadData(FolderPath,Index);
% HDW-Future
FolderPath = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\HDW\2100';
ReadData(FolderPath,Index);
% CFWI-Current
FolderPath = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\CFWI\2010';
Index = 'cfwi';
ReadData(FolderPath, Index);
% CFWI-Future
FolderPath = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\CFWI\2100';
Index = 'cfwi';
ReadData(FolderPath, Index);
% ffwi

%% Combine the 10-year time series into one matrix: each row is one pixel, each column is one date (YYYY/MM/DD).

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\HDW\2010\Annual
% cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\CFWI\2100\Annual

% dircfwi = dir('cfwi*.mat');
dirHDW = dir('HDW*.mat');
%  for i = 1:length(dircfwi)
for i = 1:length(dirHDW)
%      File = dircfwi(i).name;% load variables named "Index_year"
     File = dirHDW(i).name;% load variables named "Index_year"
     load(File);
     [r,c,h] = size(Index_year);
     npixel = r*c; % this will be the same for all files, but h may be different
     
     TempSeries = reshape(Index_year,npixel,h);
     if (i ==1 )
         TimeSeries = TempSeries;
     else
         TimeSeries = [TimeSeries,TempSeries];
     end
 end
 
 save('hdw-current.mat','TimeSeries','-v7.3');
 
 
 cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\HDW\2100\Annual
% cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\CFWI\2100\Annual

% dircfwi = dir('cfwi*.mat');
dirHDW = dir('HDW*.mat');
%  for i = 1:length(dircfwi)
for i = 1:length(dirHDW)
%      File = dircfwi(i).name;% load variables named "Index_year"
     File = dirHDW(i).name;% load variables named "Index_year"
     load(File);
     [r,c,h] = size(Index_year);
     npixel = r*c; % this will be the same for all files, but h may be different
     
     TempSeries = reshape(Index_year,npixel,h);
     if (i ==1 )
         TimeSeries = TempSeries;
     else
         TimeSeries = [TimeSeries,TempSeries];
     end
 end
 
 save('hdw-future.mat','TimeSeries','-v7.3');
 
%% plot 10-year mean to check the dataset
load C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\UnitedStates_Polygon.mat
load C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2010.mat
load C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2100.mat
load C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\GridLat_Lon;

[r,c] = size(LandFuture);
npixel = r*c;
Land = reshape(LandCurrent,npixel,1);

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices
dirInfo = dir('*.mat');
for i = 1:length(dirInfo)
    load(dirInfo(i).name);
    TimeSeries (Land==0,:) = nan;
    Avg10(:,i) = nanmean(TimeSeries,2);
    Avg = reshape(Avg10(:,i),r,c);
    subplot(3,2,i);
    h = pcolor(Lon,Lat,Avg);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');
    Variable = dirInfo(i).name(1: end-4);
    colormap(gca, jet);
    colorbar;
    xlim([-130,-65]);
    ylim([25,50]);
    title(Variable);
end

%% pre-normalization: maximum/minmum/range

load C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\UnitedStates_Polygon.mat
load C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2010.mat
load C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2100.mat
load C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\GridLat_Lon;

[r,c] = size(Lat);
npixel = r*c;

Mask2100 = reshape(LandFuture,npixel,1);
Mask2010 = reshape(LandCurrent,npixel,1);

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices
dirInfo = dir('*.mat');

load('cfwi-current.mat');
cfwi_current = TimeSeries;
load('cfwi-future.mat');
cfwi_future = TimeSeries;

AvgMaxPcnt1 = nan(npixel,length(dirInfo));
AvgMaxPcnt5 = nan(npixel,length(dirInfo));
AvgMinPcnt1 = nan(npixel,length(dirInfo));

for i = 1:length(dirInfo)
    % clean the data
    % exclude all the values out of CONUS
    load (dirInfo(i).name);
    CleanTS = TimeSeries;
    if contains(dirInfo(i).name,"current.mat")
        CleanTS(Mask2010==0,:) = nan;
    elseif contains(dirInfo(i).name,"future.mat")
        CleanTS(Mask2100==0,:) = nan;
    end
    
    % exclude all the values out of the cfwi fire season
    if(dirInfo(i).name~="cfwi*.mat")
        if contains(dirInfo(i).name,"current.mat")
             CleanTS (isnan(cfwi_current)) = nan;
        elseif contains(dirInfo(i).name,"future.mat")
             CleanTS (isnan(cfwi_future)) = nan;
        end
    end
    
    % set all the values smaller than 1 as NaN
     CleanTS ( CleanTS <1)=nan;
     
     save(['Clean-',dirInfo(i).name(1:end-4)],'CleanTS','-v7.3');
     
    % find the 1% maximum and 5% maximum value
    % exclude all the nan values (out of United States), and dates out of cfwi fire season
    CleanTS(isnan(CleanTS))=-1;
    [DescendTS,~] = sort(CleanTS,2,'descend');

    % 1% maximum
    n_1percentile = fix(0.01*size(CleanTS,2));
    percentile1_Max = DescendTS(:,1:n_1percentile);
%     DIndex_Max = DIndex(:,1:n_1percentile);
    percentile1_Max(percentile1_Max==-1) = nan;
    AvgMaxPcnt1(:,i) = nanmean(percentile1_Max,2);

    % 5% maximum
    n_5percentile = fix(0.05*size(CleanTS,2));
    percentile5_Max = DescendTS(:,1:n_5percentile);
    percentile5_Max(percentile5_Max==-1) = nan;
    AvgMaxPcnt5(:,i) = nanmean(percentile5_Max,2);
    save( ['Max-',dirInfo(i).name(1:end-4)],'percentile5_Max','percentile1_Max','-v7.3');
    clear DescendTS percentile1_Max percentile5_Max
    
    % find current min for normalization
    if contains(dirInfo(i).name,"current.mat")
        % find the 1% minimum
        CleanTS(CleanTS==-1)=9999;
        [AscendTS,~] = sort(CleanTS,2,'ascend');
        percentile1_Min = AscendTS(:,1:n_1percentile);
        percentile1_Min(percentile1_Min==9999) = nan;
        AvgMinPcnt1(:,i) = nanmean(percentile1_Min,2);
        save( ['Min-',dirInfo(i).name(1:end-4)],'percentile1_Min','-v7.3');
        clear AscendTS percentile1_Min
    end
     
end

% calculate the range 
% AvgMaxPcnt1(AvgMaxPcnt1==-1) = NaN;
% AvgMaxPcnt5(AvgMaxPcnt5==-1) = NaN;
% AvgMinPcnt1(AvgMinPcnt1==9999) = NaN;
Range = AvgMaxPcnt1-AvgMinPcnt1;
Range = [Range(:,1),Range(:,3)];

save('MaxPcnt1.mat','AvgMaxPcnt1','-v7.3');
save('MaxPcnt5.mat','AvgMaxPcnt5','-v7.3');
save('Range.mat','Range','-v7.3');

%% prepare files for normalization
clear

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices\Analysis

load('Range.mat');
load('Min-cfwi-current.mat');
load('Max-cfwi-current.mat');
AvgMinPcnt1 = nanmean(percentile1_Min,2);
AvgMaxPcnt1 = nanmean(percentile1_Max,2);
AvgMaxPcnt1(AvgMaxPcnt1==-1) = NaN;
AvgMinPcnt1(AvgMinPcnt1==9999) = NaN;
range = Range(:,1);
save('cfwi-stat.mat','AvgMinPcnt1','AvgMaxPcnt1','range','-v7.3');

clear
load('Range.mat');
load('Min-hdw-current.mat');
load('Max-hdw-current.mat');
AvgMinPcnt1 = nanmean(percentile1_Min,2);
AvgMaxPcnt1 = nanmean(percentile1_Max,2);
AvgMaxPcnt1(AvgMaxPcnt1==-1) = NaN;
AvgMinPcnt1(AvgMinPcnt1==9999) = NaN;
range = Range(:,2);
save('hdw-stat.mat','AvgMaxPcnt1','AvgMinPcnt1','range','-v7.3');

%% Normalization
clear

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices\Analysis
% load Clean Time Series: exclude all pixels out of CONUS/fire season/smaller than 1 
dirInfo = dir('Clean*.mat');

for i = 3:length(dirInfo)
    
    load(dirInfo(i).name);
    
    if contains(dirInfo(i).name, "hdw")
       load('hdw-stat.mat');
       Std_TS = (CleanTS-AvgMinPcnt1)./range;
       
    elseif contains(dirInfo(i).name,"cfwi")
       load('cfwi-stat.mat');
       Std_TS = (CleanTS-AvgMinPcnt1)./range; 
    end 
    
     save( ['N-',dirInfo(i).name(7:end-4)],'CleanTS','Std_TS','-v7.3');
end

%% count number/percentage of dates when normalized value >0.5 bad days

clear 

load C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\GridLat_Lon;
[r,c] = size(Lat);
npixel = r*c;

% number of extreme/bad days
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices\Analysis
dirInfo = dir('N-*.mat');
Bad = nan(npixel,length(dirInfo));
Extreme = nan(npixel,length(dirInfo));

for i = 1:length(dirInfo)
    disp(i);
    load(dirInfo(i).name);
    T(i) = string(dirInfo(i).name(3:end-4));
    BadDays= Std_TS>0.5;
    Bad(:,i) = nansum(BadDays,2);
    ExtremeDays = Std_TS>0.95;
    Extreme(:,i) = nansum(ExtremeDays,2);
end

Extreme(Extreme==0) = nan;
Bad(Bad==0) = nan;

% length of the cfwi fire season
load ('N-cfwi-future.mat');
[r,c] = size(Lat);
Index = CleanTS;
Index (~isnan(Index)) = 1;
cfwiFuture_FL = nansum(Index,2);

load ('N-cfwi-current.mat');
Index = CleanTS;
Index (~isnan(Index)) = 1;
cfwiCurrent_FL = nansum(Index,2);

cfwi_FL = [cfwiCurrent_FL,cfwiFuture_FL];
cfwi_FL(cfwi_FL==0) = nan;

%% count the number of consecutive bad days
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices\Analysis
dirInfo = dir('N-*.mat');

for i = 3:length(dirInfo)
    
    disp(i);
    load(dirInfo(i).name);
    Index = (Std_TS>0.5)';
    fiDiff = diff(Index);

    % check whether this scenario exists: the first day of the time series >0.5
    a = sum(Index(1,:)>0.5);% if a = 0, then no such case
    disp(a);
    
    % cfwi doesn't have this case
    fiDiff = [zeros(1,length(fiDiff));fiDiff];

    % Idea: for each pixel and each date in matrix cfwiDiff (each column is 
    % one pixel, each row is one date), calculate the sum of its value and
    % its following four days. If the sum is 1 and the value of that date is 1,
    % then it's consecutive bad days

    % matrix to calculate the sum of five days
    Transform = zeros(size(fiDiff,1));
    n = size(Transform,1);
    for j = 1:n

        if j >=2 && j<(n-4)
            Transform(j,j:(j+4))=1;
        end
    end


    % Test as an example
    % a = Transform(1:8,1:8);
    % 
    % b1 = ones(1,10);
    % b2 = 2.*b1;
    % b3 = 3.*b1;
    % b4 = 4.*b1;
    % b5 = 5.*b1;
    % b6 = 6.*b1;
    % b7 = 7.*b1;
    % b8 = 8.*b1;
    % 
    % b = [b1;b2;b3;b4;b5;b6;b7;b8];
    % c =a*b;

    % sum the 5-day difference
    fiSum5 = Transform*fiDiff;


    % a = [cfwiDiff(1341:1359,:);cfwiSum5(1341:1359,:)];

    Value = (fiDiff==1)';
    Sum5 = (fiSum5==1)';
    %  a = [cfwiDiff(1341,:);cfwiSum5(1341,:);Value(1341,:);Sum(1341,:)];
    Cnsct = Value==1 & Sum5==1;
    nCnsct = sum(Cnsct,2);
    nCnsct(nCnsct==0) = nan;

    save(['Cnsct-',dirInfo(i).name(1:end-4)],'fiDiff','fiSum5','nCnsct','-v7.3');
end

a= reshape(nCnsct,r,c);
h = pcolor(Lon,Lat,a);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

%% Student t-test for Max 5%
% cfwi
load('Max-cfwi-future.mat');
futurecfwi = percentile5_Max;

load('Max-cfwi-current.mat');
currentcfwi = percentile5_Max;

clear percentile1_Max percentile5_Max

h_cfwi = nan(size(currentcfwi,1),1);
p_cfwi = nan(size(currentcfwi,1),1);

for i = 1:length(currentcfwi)
    [h_cfwi(i,1),p_cfwi(i,1)] = ttest2(futurecfwi(i,:),currentcfwi(i,:));
end

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices\Analysis
% HDW
load('Max-hdw-future.mat');
futurehdw = percentile5_Max;

load('Max-hdw-current.mat');
currenthdw = percentile5_Max;

clear percentile1_Max percentile5_Max

h_hdw = nan(size(currenthdw,1),1);
p_hdw = nan(size(currenthdw,1),1);

for i = 1:length(currenthdw)
    [h_hdw(i,1),p_hdw(i,1)] = ttest2(futurehdw(i,:),currenthdw(i,:));
end

save('ttest.mat','h_cfwi','p_cfwi','h_hdw','p_hdw');
