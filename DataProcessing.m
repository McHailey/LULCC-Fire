%%% This is a script that get/read/analyze WRF outputs
%% read .nc files into matrix and group by year/month

% data under future LULC
Folder = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2100\Raw';
ReadData(Folder);

% data under current LULC
Folder = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2010\Raw';
ReadData(Folder);
%% Calculate average annual cycle (daily)
% current
AvgCycleOutput = 'AvgCycle2010';
AnnualPath = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2010\Annual';
Avgffwi_Current = AvgCycle(AnnualPath,AvgCycleOutput);

% future
AvgCycleOutput = 'AvgCycle2100';
AnnualPath = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2100\Annual';
Avgffwi_Future = AvgCycle(AnnualPath,AvgCycleOutput);

%% calculate the average annual cycle for current land cover (seasonally)

% Dec-Feb: Winter; March-May: Spring; June-August: Summer;
% September-November: Autumn
% average in season

% % current
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2010\Monthly
% 
% % Winter
% CurrentWinterDir = [dir('ffwi*12.mat'); dir('ffwi*01.mat'); dir('ffwi*02.mat')];
% % % Spring
% CurrentSpringDir = [dir('ffwi*03.mat'); dir('ffwi*04.mat'); dir('ffwi*05.mat')];
% % % Summer
% CurrentSummerDir = [dir('ffwi*06.mat'); dir('ffwi*07.mat'); dir('ffwi*08.mat')];
% % % Autumn
% CurrentAutumnDir = [dir('ffwi*09.mat'); dir('ffwi*10.mat'); dir('ffwi*11.mat')];
% CurrentWinterAvg2D = TemporalAvg(CurrentWinterDir,'WinterAvg2D');
% CurrentSpringAvg2D = TemporalAvg(CurrentSpringDir,'SpringAvg2D');
% CurrentSummerAvg2D = TemporalAvg(CurrentSummerDir,'SummerAvg2D');
% CurrentAutumnAvg2D = TemporalAvg(CurrentAutumnDir,'AutumnAvg2D');

load('WinterAvg2D.mat');
CurrentWinter = Ave2D;

load('SpringAvg2D.mat');
CurrentSpring = Ave2D;

load('SummerAvg2D.mat');
CurrentSummer = Ave2D;

load('AutumnAvg2D.mat');
CurrentAutumn= Ave2D;

CurrentSeason = cat(3,CurrentSpring,CurrentSummer,CurrentAutumn,CurrentWinter);

%%%%% -----------------------------Future-----------------------------%%%%%
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2100\Monthly;
% Winter
% FutureWinterDir = [dir('ffwi*12.mat'); dir('ffwi*01.mat'); dir('ffwi*02.mat')];
% % Spring
% FutureSpringDir = [dir('ffwi*03.mat'); dir('ffwi*04.mat'); dir('ffwi*05.mat')];
% % Summer
% FutureSummerDir = [dir('ffwi*06.mat'); dir('ffwi*07.mat'); dir('ffwi*08.mat')];
% % Autumn
% FutureAutumnDir = [dir('ffwi*09.mat'); dir('ffwi*10.mat'); dir('ffwi*11.mat')];
% FutureWinterAvg2D = TemporalAvg(FutureWinterDir,'WinterAvg2D');
% FutureSpringAvg2D = TemporalAvg(FutureSpringDir,'SpringAvg2D');
% FutureSummerAvg2D = TemporalAvg(FutureSummerDir,'SummerAvg2D');
% FutureAutumnAvg2D = TemporalAvg(FutureAutumnDir,'AutumnAvg2D');

load('WinterAvg2D.mat');
FutureWinter = Ave2D;

load('SpringAvg2D.mat');
FutureSpring = Ave2D;

load('SummerAvg2D.mat');
FutureSummer = Ave2D;

load('AutumnAvg2D.mat');
FutureAutumn= Ave2D;

FutureSeason = cat(3,FutureSpring,FutureSummer,FutureAutumn,FutureWinter);

%%%% ----------------------------------Masking------------------------%%%%
% load Lat/Lon for data grids
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project
load('GridLat_Lon.mat');
% load mask indices
% load('MaskUS.mat');
load('Mask2010.mat');
load('Mask2100.mat');

%Current

CleanCurrentSpring = MaskPoly(Lat,Lon,CurrentSpring,LandCurrent);
CleanCurrentSummer = MaskPoly(Lat,Lon,CurrentSummer,LandCurrent);
CleanCurrentAutumn = MaskPoly(Lat,Lon,CurrentAutumn,LandCurrent);
CleanCurrentWinter = MaskPoly(Lat,Lon,CurrentWinter,LandCurrent);

%Future
CleanFutureSpring = MaskPoly(Lat,Lon,FutureSpring,LandFuture);
CleanFutureSummer = MaskPoly(Lat,Lon,FutureSummer,LandFuture);
CleanFutureAutumn = MaskPoly(Lat,Lon,FutureAutumn,LandFuture);
CleanFutureWinter = MaskPoly(Lat,Lon,FutureWinter,LandFuture);

% Difference
SpringDifference = CleanFutureSpring - CleanCurrentSpring;
SummerDifference = CleanFutureSummer - CleanCurrentSummer;
AutumnDifference = CleanFutureAutumn - CleanCurrentAutumn;
WinterDifference = CleanFutureWinter - CleanCurrentWinter;

%% calculate the average value with current/future LULC

CleanCurrent = cat(3,CleanCurrentSpring,CleanCurrentSummer, CleanCurrentAutumn, CleanCurrentWinter);
MeanCurrent = mean(CleanCurrent,3);

CleanFuture = cat(3,CleanFutureSpring,CleanFutureSummer, CleanFutureAutumn, CleanFutureWinter);
MeanFuture = mean(CleanFuture,3); 

[y,x] = size(MeanFuture);

% write annual average surface to nc file
nccreate('ffwi_Current.nc','AnnualMean','Dimensions', {'y',y,'x',x});
ncwrite('ffwi_Current.nc','AnnualMean',MeanCurrent);
nccreate('ffwi_Future.nc','AnnualMean','Dimensions', {'y',y,'x',x});
ncwrite('ffwi_Future.nc','AnnualMean',MeanFuture);


h = pcolor(Lon, Lat, MeanCurrent);
h.EdgeColor = 'none';
h.ZData = h.CData;
title('ffwi 10 year Average (Current)');
xlim ([-130,-65]);
ylim([23,52]);colorbar;


%%%% -------------------------Seasonal Plotting----------------------- %%%%
%Current
figure;
% mstruct = defaultm('lambert');
% mstruct.geoid = wgs84Ellipsoid;
% mstruct = defaultm(mstruct);
% 
% [y_V,x_V] =mfwdtran(mstruct,Lat_V,Lon_V);
% y = reshape(y_V,n,m);
% x = reshape(x_V,n,m);

% Current Spring
subplot(4,3,1);
contourf(Lon,Lat,CleanCurrentSpring,'LineStyle','none');
% h1 = pcolor(Lon,Lat,CleanCurrentSpring);
% h1.EdgeColor = 'none';
% h1.ZData = h1.CData;
title('Spring 10 year Average (Current)');
c1 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([10,60]);
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

% Future Spring
subplot(4,3,2);
% h1 = pcolor(Lon,Lat,CleanFutureSpring);
% h1.EdgeColor = 'none';
% h1.ZData = h1.CData;
contourf(Lon,Lat,CleanFutureSpring,'LineStyle','none');
title('Spring 10 year Average (Future)');
c5 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([10,60]);
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

% Differences in Spring
subplot(4,3,3);
% h1 = pcolor(Lon,Lat,SpringDifference);
contourf(Lon,Lat,SpringDifference,'LineStyle','none');
% h1.EdgeColor = 'none';
% h1.ZData = h1.CData;
title('Differencesin Spring(Future-Current)');
c9 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([-10,5]);
colormap('jet');
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

% Current Summer
subplot(4,3,4);
% figure;
contourf(Lon,Lat,CleanCurrentSummer,'LineStyle','none');
title('Summer 10 year Average (Current)');
c2 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([10,60]);
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');



% Future Summer
subplot(4,3,5);
contourf(Lon,Lat,CleanFutureSummer,'LineStyle','none');
% h2 = pcolor(Lon,Lat,CleanFutureSummer);
% h2.EdgeColor = 'none';
% h2.ZData = h2.CData;
title('Summer 10 year Average (Future)');
c6 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([10,60]);
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

subplot(4,3,6);
contourf(Lon,Lat,SummerDifference,'LineStyle','none');
% h1 = pcolor(Lon,Lat,SummerDifference);
% h1.EdgeColor = 'none';
% h1.ZData = h1.CData;
title('Differences in Summer(Future-Current)');
c10 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([-10,5]);
colormap('jet');
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

% Current Autumn
subplot(4,3,7);
contourf(Lon,Lat,CleanCurrentAutumn,'LineStyle','none');
% h3 = pcolor(Lon,Lat,CleanCurrentAutumn);
% h3.EdgeColor = 'none';
% h3.ZData = h3.CData;
title('Autumn 10 year Average (Current)');
c3 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([10,60]);
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

% Future Autumn
subplot(4,3,8);
contourf(Lon,Lat,CleanFutureAutumn,'LineStyle','none');
% h3 = pcolor(Lon,Lat,CleanFutureAutumn);
% h3.EdgeColor = 'none';
% h3.ZData = h3.CData;
title('Autumn 10 year Average (Future)');
c7 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([10,60]);
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

subplot(4,3,9);
contourf(Lon,Lat,AutumnDifference,'LineStyle','none');
% h1 = pcolor(Lon,Lat,AutumnDifference);
% h1.EdgeColor = 'none';
% h1.ZData = h1.CData;
title('Differences in Autumn(Future-Current)');
c11 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([-10,5]);
colormap('jet');
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

% Current Winter
subplot(4,3,10);
contourf(Lon,Lat,CleanCurrentWinter,'LineStyle','none');
% h4 = pcolor(Lon,Lat,CleanCurrentWinter);
% h4.EdgeColor = 'none';
% h4.ZData = h4.CData;
title('Winter 10 year Average (Current)');
c4 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([10,60]);
% colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

% Future Winter
subplot(4,3,11);
contourf(Lon,Lat,CleanFutureWinter,'LineStyle','none');
% h4 = pcolor(Lon,Lat,CleanFutureWinter);
% h4.EdgeColor = 'none';
% h4.ZData = h4.CData;
title('Winter 10 year Average (Future)');
c8 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([10,60]);
colorbar('Position',[0.08 0.10 0.015 0.8],'LineWidth',1);
hold on;
plot([states.Lon],[states.Lat], 'Color','black');

subplot(4,3,12);
contourf(Lon,Lat,WinterDifference,'LineStyle','none');
% h1 = pcolor(Lon,Lat,WinterDifference);
% h1.EdgeColor = 'none';
% h1.ZData = h1.CData;
title('Differences in Winter(Future-Current)');
c12 = caxis;
xlim ([-130,-65]);
ylim([25,50]);
caxis([-10,5]);
colormap('jet');
colorbar;
hold on;
plot([states.Lon],[states.Lat], 'Color','black');


c1 = [min([c1,c2,c3,c4,c5,c6,c7,c8]),max([c1,c2,c3,c4,c5,c6,c7,c8])];
colorbar;
c2= [min([c9,c10,c11,c12]),max([c9,c10,c11,c12])];
colorbar;

% Seasonal Boxplot (Average 10 year, Average region)
CurrentWinterAvg1D = reshape(CleanCurrentWinter,92981,1);
CurrentSpringAvg1D = reshape(CleanCurrentSpring,92981,1);
CurrentAutumnAvg1D = reshape(CleanCurrentAutumn,92981,1);
CurrentSummerAvg1D = reshape(CleanCurrentSummer,92981,1);
CurrentSeasonalVector = [CurrentSpringAvg1D, CurrentSummerAvg1D,CurrentAutumnAvg1D,CurrentWinterAvg1D];

figure;

subplot(1,2,1);
boxplot(CurrentSeasonalVector);
set(gca,'XTickLabel',{'MAM(Spring)','JJA(Summer)','SON(Autumn)','DJF(Winter)'});
title('Seasonal Boxplot - 10 year with Current LULC');

% Seasonal Boxplot (Average 10 year, Average region)
FutureWinterAvg1D = reshape(CleanFutureWinter,92981,1);
FutureSpringAvg1D = reshape(CleanFutureSpring,92981,1);
FutureAutumnAvg1D = reshape(CleanFutureAutumn,92981,1);
FutureSummerAvg1D = reshape(CleanFutureSummer,92981,1);
FutureSeasonalVector = [FutureSpringAvg1D, FutureSummerAvg1D,FutureAutumnAvg1D,FutureWinterAvg1D];

subplot(1,2,2);
boxplot(FutureSeasonalVector);
set(gca,'XTickLabel',{'MAM(Spring)','JJA(Summer)','SON(Autumn)','DJF(Winter)'});
title('Seasonal Boxplot - 10 year with Future LULC');

figure
histogram(CurrentSeasonalVector,'facealpha',.5,'BinWidth',1);
title('Histogram for ffwi in Current LULC');

hold on;

figure
histogram(FutureSeasonalVector,'facealpha',.5,'BinWidth',1);
title('Histogram for ffwi in Future LULC');