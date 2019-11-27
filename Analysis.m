%%% This script is used to do preliminary analysis for Forsberg Fire indices
%%% with current LULC map. 

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2010
%% calculate the average 10-year annual value for current land cover
% load('ffwi_current.mat');
% ffwi_AnnualMN = mean(ffwi,3);

%% calculate the average annual cycle for current land cover (daily)

% for 2/29, only leap years have it
x = [2019 2 29];
x_date = datetime(x);
x_doy = day(x_date,'dayofyear');

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2010\Annual
dirInfo = dir('*.mat');
j = 0;
for i = 1:length(dirInfo)
    FileName = dirInfo(i).name;
    load(FileName);    
    if (mod (str2num(FileName(end-7:end-4)),4)==0)
        
        % save data for 02/29 of each year
        if j==0
          Leap = ffwi(:,:,x_doy);
        else
          Leap = cat(3,Leap,ffwi(:,:,x_doy));
        end
        % clean data for 02/29 in the annually stacked dataset
        ffwi(:,:,x_doy) = [];
        j = j+1;
    end
    
        % each file has a 259*359*364 matrix (daily stack)
    [l,w,h] = size(ffwi);
    % annual stack: each year is a layer
    % each column is one day for the whole region
    % each row is the annual cycle of one pixel
    ffwi_v(:,:,i) = reshape(ffwi,l*w,h);
end

% calculate the average daily annual cycle (without 02/29)
mean(ffwi_v,3);
save('ffwi_v')

%% calculate the average annual cycle for current land cover (seasonally)
% Dec-Feb: Winter; March-May: Spring; June-August: Summer;
% September-November: Autumn
% average in season
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2010\Monthly
% Winter
WinterDir = [dir('ffwi*12.mat'); dir('ffwi*01.mat'); dir('ffwi*02.mat')];
% Spring
SpringDir = [dir('ffwi*03.mat'); dir('ffwi*04.mat'); dir('ffwi*05.mat')];
% Summer
SummerDir = [dir('ffwi*06.mat'); dir('ffwi*07.mat'); dir('ffwi*08.mat')];
% Autumn
AutumnDir = [dir('ffwi*09.mat'); dir('ffwi*10.mat'); dir('ffwi*11.mat')];
WinterAvg2D = TemporalAvg(WinterDir,'WinterAvg2D');
SpringAvg2D = TemporalAvg(SpringDir,'SpringAvg2D');
SummerAvg2D = TemporalAvg(SummerDir,'SummerAvg2D');
AutumnAvg2D = TemporalAvg(AutumnDir,'AutumnAvg2D');


subplot(2,2,1);
h1 = imagesc(SpringAvg2D);
h1.EdgeColor = 'none';
h1.ZData = h1.CData;
title('Spring 10 year Average (Current)');

subplot(2,2,2);
h2 = pcolor(SummerAvg2D);
h2.EdgeColor = 'none';
h2.ZData = h2.CData;
title('Summer 10 year Average (Current)');

subplot(2,2,3);
h3 = pcolor(AutumnAvg2D);
h3.EdgeColor = 'none';
h3.ZData = h3.CData;
title('Autumn 10 year Average (Current)');

subplot(2,2,4);
h4 = pcolor(WinterAvg2D);
h4.EdgeColor = 'none';
h4.ZData = h4.CData;
title('Winter 10 year Average (Current)');

% Seasonal Boxplot (Average 10 year, Average region)
WinterAvg1D = reshape(WinterAvg2D,92981,1);
SpringAvg1D = reshape(SpringAvg2D,92981,1);
AutumnAvg1D = reshape(AutumnAvg2D,92981,1);
SummerAvg1D = reshape(SummerAvg2D,92981,1);
SeasonalVector = [SpringAvg1D, SummerAvg1D,AutumnAvg1D,WinterAvg1D];
boxplot(SeasonalVector);
set(gca,'XTickLabel',{'MAM(Spring)','JJA(Summer)','SON(Autumn)','DJF(Winter)'});

%% calculate the average value with current/future LULC

CleanCurrent = cat(3,CleanCurrentSpring,CleanCurrentSummer, CleanCurrentAutumn, CleanCurrentWinter);
MeanCurrent = mean(CleanCurrent,3);

CleanFuture = cat(3,CleanFutureSpring,CleanFutureSummer, CleanFutureAutumn, CleanFutureWinter);
MeanFuture = mean(CleanFuture,3); 

% write annual average surface to nc file
nccreate('ffwi_Current.nc','AnnualMean','Dimensions', {'x',l,'y',w});
ncwrite('ffwi_Current.nc','AnnualMean',MeanCurrent);

h = pcolor(Lon, Lat, MeanCurrent);
h.EdgeColor = 'none';
h.ZData = h.CData;
title('ffwi 10 year Average (Current)');
xlim ([-130,-65]);
ylim([23,52]);colorbar;

nccreate('ffwi_Future.nc','AnnualMean','Dimensions', {'x',l,'y',w});
ncwrite('ffwi_Future.nc','AnnualMean',MeanFuture);
%% calculate the average value for a specfic year/month

load('ffwi200701.mat');
[l,w,h] = size(ffwi_month);
npixel = l*w;
ffwi2007_01 = reshape(ffwi_month,npixel,h);
MN_ffwi200701 = mean(ffwi2007_01,2);
MN_ffwi = reshape(MN_ffwi200701,l,w);
save('MN_ffwi200701','MN_ffwi');
% write nc files
nccreate('ffwi200601.nc','Mean','Dimensions', {'y',l,'x',w});
ncwrite('ffwi200601.nc','Mean',MN_ffwi_200601);
nccreate('ffwi200701.nc','Mean','Dimensions', {'y',l,'x',w});
ncwrite('ffwi200701.nc','Mean',MN_ffwi_200701);

%% Generate a mask that only includes land masses in CONUS
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project
load('LT2010.mat');
load('LT2100.mat');

% Mask of current land cover type
LandCurrent = LT2010_Grid~=1 & MaskIndex==1;
% Mask of future land cover type
LandFuture = LT2100_Grid~=1 & MaskIndex==1;

