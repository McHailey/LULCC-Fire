%  cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2100\Annual
% Folder = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\HDW\2100\Annual';
% Folder = 'C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\CFWI\2100\Annual';
cd(Folder);
 %% Combine the 10-year time series into one matrix: each row is one pixel, each column is one date (YYYY/MM/DD).
 dircfwi = dir('cfwi*.mat');
 for i = 1:length(dircfwi)
     File = dircfwi(i).name;% load variables named "Index_year"
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
 
% save('10YearTS.mat','TimeSeries');
% cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\HDW2100\Annual
% load ('10YearTS.mat');
% load('ffwi_FutureTS.mat');
%% sort each row in the matrix
Clean_TS = TimeSeries;
% for cfwi, On the dates out of the growing season, the pixel will be set as NaN.
% Change it to -1 for sorting the row data. 
Clean_TS(isnan(Clean_TS))=-1;
[Descend_TS,~] = sort(Clean_TS,2,'descend');

%% find the extreme maximum value for each pixel (each row)
% Current
%  load('C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2010.mat');

%Future
load('C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2100.mat');
[r,c] = size(LandFuture);
npixel = r*c;

% Extreme_Max = reshape(Descend_TS(:,1),r,c);
% Extreme_Max(LandFuture==0) = nan;

%% find the 5% maximum values for each pixel (each row)
n_5percentile = fix(0.05*size(TimeSeries,2));
percentile5_Max = Descend_TS(:,1:n_5percentile);
mean_percentile5 = mean(percentile5_Max,2);
Pcnt5_Max = reshape( mean_percentile5,r(1),c(1));
Pcnt5_Max(LandFuture==0) = nan;

%% find the 1% maximum 
n_1percentile = fix(0.01*size(TimeSeries,2));
percentile1_Max = Descend_TS(:,1:n_1percentile);
mean_percentile1 = nanmean(percentile1_Max,2);
Pcnt1_Max = reshape( mean_percentile1,r,c);
% Pcnt1_Max (LandCurrent==0) = nan;
Pcnt1_Max (LandFuture==0) = nan;
clear Descend_TS TempSeries
% %% find the 10% maximum
% n_10percentile = fix(0.1*size(TimeSeries,2));
% percentile10_Max = Descend_TS(:,1:n_10percentile);
% mean_percentile10 = mean(percentile10_Max,2);
% Pcnt10_Max = reshape( mean_percentile10,r(1),c(1));
% Pcnt10_Max (LandFuture==0) = nan;

% %% make plots (maximum values)
% load('UnitedStates_Polygon.mat');
% load('GridLat_Lon.mat');
% 
% figure;
% subplot(2,2,1);
% h1 = pcolor(Lon,Lat, Extreme_Max);
% h1.EdgeColor = 'none';
% % h1.ZData =h1.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% title('Maximum ffwi');
% % caxis([30,120]);
% colorbar;
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');


% subplot(2,2,2)
% h2 = pcolor(Lon,Lat, Pcnt1_Max);
% h2.EdgeColor = 'none';
% % h2.ZData =h2.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% title('highest 1% mean, ffwi');
% % caxis([30,120]);
% colorbar;
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% 
% subplot(2,2,3);
% h3 = pcolor(Lon,Lat, Pcnt5_Max);
% h3.EdgeColor = 'none';
% % h3.ZData =h3.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% title('highest 5% mean');
% % caxis([30,120]);
% colorbar;
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');


% subplot(2,2,4)
% h4 = pcolor(Lon,Lat, Pcnt10_Max);
% h4.EdgeColor = 'none';
% % h4.ZData =h4.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% title('highest 10% mean');
% % caxis([30,120]);
% colorbar;

%% sort each row in the matrix
% for cfwi, to improve computational efficiency, set -1 (previously
% converted from NaN) to 9999;

Clean_TS(Clean_TS<1) = 9999;% Values<0.1 doesn't matter any way

[Ascend_TS,~] = sort(Clean_TS,2,'ascend');

%% find the extreme minimum value for each pixel (each row)
% load('Mask2010.mat');
% Extreme_Min = reshape(Ascend_TS(:,1),r(1),c(1));
% Extreme_Min(LandFuture==0) = nan;

% %% find the 5% minimum values for each pixel (each row)
% n_5percentile = fix(0.05*size(TimeSeries,2));
% percentile5_Min = Ascend_TS(:,1:n_5percentile);
% mean_percentile5 = mean(percentile5_Min,2);
% Pcnt5_Min = reshape( mean_percentile5,r(1),c(1));
% Pcnt5_Min(LandFuture==0) = nan;

%% find the 1% minimum values for each pixel (each row)

n_1percentile = fix(0.01*size(TimeSeries,2));
percentile1_Min = Ascend_TS(:,1:n_1percentile);
mean_percentile1 = nanmean(percentile1_Min,2);
Pcnt1_Min = reshape( mean_percentile1,r,c);
% Pcnt1_Min (LandCurrent==0) = nan;
Pcnt1_Min (LandFuture==0) = nan;

clear Ascend_TS Index_year Index_A Index_D
% %% find the 10% minimum values for each pixel (each row)
% n_10percentile = fix(0.1*size(TimeSeries,2));
% percentile10_Min = Ascend_TS(:,1:n_10percentile);
% mean_percentile10 = mean(percentile10_Min,2);
% Pcnt10_Min = reshape( mean_percentile10,r(1),c(1));
% Pcnt10_Min (LandFuture==0) = nan;

% %% plot the minimum map
% figure;
% subplot(2,2,1);
% h1 = pcolor(Lon,Lat, Extreme_Min);
% h1.EdgeColor = 'none';
% % h1.ZData =h1.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% title('Minimum ffwi');
% % caxis([30,120]);
% colorbar;
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');

% subplot(2,2,2)
% h2 = pcolor(Lon,Lat, Pcnt1_Min);
% h2.EdgeColor = 'none';
% % h2.ZData =h2.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% title('lowest 1% mean');
% % caxis([30,120]);
% colorbar;
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');

% subplot(2,2,3);
% h3 = pcolor(Lon,Lat, Pcnt5_Min);
% h3.EdgeColor = 'none';
% % h3.ZData =h3.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% title('lowest 5% mean');
% % caxis([30,120]);
% colorbar;
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');

% subplot(2,2,4)
% h4 = pcolor(Lon,Lat, Pcnt10_Min);
% h4.EdgeColor = 'none';
% % h4.ZData =h4.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% title('lowest 10% mean');
% % caxis([30,120]);
% colorbar;
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');

%% calculate the range

% exclude NaN/-1/9999 values
percentile1_Max(percentile1_Max==-1) = NaN;
percentile1_Min(percentile1_Min==9999) = NaN;
range = Pcnt1_Max-Pcnt1_Min;

%% normalization
range_V = reshape(range, npixel,1);
Pcnt1_MinV = reshape(Pcnt1_Min,npixel,1);
TimeSeries_Stdize = (TimeSeries-Pcnt1_MinV)./range_V;
Avg10 = reshape (nanmean(TimeSeries_Stdize,2),259,359);

%% Plots with normalized indices (3*2 figure plot)
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Normalization;
load('C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\UnitedStates_Polygon.mat');
dirInfo = dir('*.mat');
n = length(dirInfo);
npixel = r*c;
% figure;
T = ["HDW-Current", "HDW-Future","cfwi-Current", "cfwi-Future","ffwi-Current", "ffwi-Future"];
for i = 1 :n
    % calculate number of days > 0.5
    load(dirInfo(i).name);
    Extreme_Index = TimeSeries_Stdize>0.5;
    cumulateExtreme(:,i) = sum (Extreme_Index,2);
    Max(:,i) = reshape(Pcnt1_Max,npixel,1);
    Min(:,i) = reshape(Pcnt1_Min,npixel,1);
    Range(:,i) = reshape(range,npixel,1);

    % plot mean
%     subplot(2,4,i);
%     h = pcolor(Lon,Lat, Avg10);
%     h.EdgeColor = 'none';
%     % h4.ZData =h4.CData; 
%     hold on;
%     plot([states.Lon],[states.Lat], 'Color','black');
%     title(strcat(T(i), ", 10year Average"));
%     caxis([0.1,0.5]);
%     colorbar;
%     hold on;
%     plot([states.Lon],[states.Lat], 'Color','black');
    cumulateExtreme(cumulateExtreme==0) = NaN;
%     ExtremeDate = reshape(cumulateExtreme(:,i),r,c);
%     subplot(3,2,i);
%     h = pcolor(Lon,Lat,ExtremeDate);
%     h.EdgeColor = 'none';
%     % h4.ZData =h4.CData; 
%     hold on;
%     plot([states.Lon],[states.Lat], 'Color','black');
%     title(T(i));
% %    caxis([0.1,0.5]);
%     colorbar;
%     hold on;
%     plot([states.Lon],[states.Lat], 'Color','black');
%     save(strcat(T(i),"ExtDay.mat"), 'ExtremeDate');

end

%% Compute the fire season length,percentage of days >0.5

load ('N_cfwi_futureC.mat');
[r,c] = size(Lat);
Index = TimeSeries_Stdize;
Index (~isnan(Index)) = 1;
cfwiFuture_Length = nansum(Index,2);

load ('N_cfwi_current.mat');
Index = TimeSeries_Stdize;
Index (~isnan(Index)) = 1;
cfwiCurrent_Length = nansum(Index,2);
cfwi_Length = [cfwiCurrent_Length,cfwiFuture_Length];% 1-Current; 2-Future

constantLength = repmat(3642*ones(size(cfwiCurrent_Length)),1,2);
FS_Length = [constantLength,cfwi_Length,constantLength];

load('NumExtremeDay.mat');
ExtPect = 100.*cumulateExtreme./FS_Length;

figure;
for i = 1:size(ExtPect,2)
    subplot(3,2,i);
    Pect = reshape(ExtPect(:,i),r,c);
    h = pcolor(Lon,Lat,Pect);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');
    colorbar;
   
end

figure;
for i = 1:size(ExtPect,2)
    subplot(3,2,i);
    Pect = reshape(Min(:,i),r,c);
    h = pcolor(Lon,Lat,Pect);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');
    colorbar;
   
end

%% If use Canadian Fire Season on HDW/FFWI

% future
load ('N_cfwi_futureC.mat');
cfwi = TimeSeries_Stdize;
load('N_HDW_futureC.mat');
hdw = TimeSeries_Stdize;
hdw(isnan(cfwi)) = nan;
% load('N_ffwi_FutureC.mat');
% ffwi = TimeSeries_Stdize;
% ffwi(isnan(cfwi)) = nan;

Extreme_hdw = hdw>0.95;
Cumulate_hdw = nansum(Extreme_hdw,2);
Extreme_cfwi = cfwi>0.95;
Cumulate_cfwi = nansum(Extreme_cfwi,2);

% 
% Bad_hdw = hdw>0.5;
% Cumulate_hdw = nansum(Bad_hdw,2);
% Bad_cfwi = cfwi>0.5;
% Cumulate_cfwi = nansum(Bad_cfwi,2);

load('NumExtremeDay.mat');
Pcnt_ffwi_Future2 = 100.*Cumulate_ffwi./cfwi_Length(:,2);% cfwi_length(:,2) is for future
Pcnt_hdw_Future2 = 100.*Cumulate_hdw./cfwi_Length(:,2);

% current

load ('N_cfwi_current.mat');
cfwi = TimeSeries_Stdize;
load('N_HDW_current.mat');
hdw = TimeSeries_Stdize;
hdw(isnan(cfwi)) = nan;
% load('N_ffwi_current.mat');
% ffwi = TimeSeries_Stdize;
% ffwi(isnan(cfwi)) = nan;

Extreme_hdwC = hdw>0.95;
Cumulate_hdwC = nansum(Extreme_hdwC,2);
Extreme_cfwiC = cfwi>0.95;
Cumulate_cfwiC = nansum(Extreme_cfwiC,2);

ExtremeNum = [Cumulate_hdwC,Cumulate_hdw,Cumulate_cfwiC,Cumulate_cfwi];


load('NumExtremeDay.mat');
Pcnt_ffwi_Current2 = 100.*Cumulate_ffwi./cfwi_Length(:,1);% cfwi_length(:,2) is for future
Pcnt_hdw_Current2 = 100.*Cumulate_hdw./cfwi_Length(:,1);

% Plot
subplot(2,2,1);
Pect = reshape(Pcnt_hdw_Current2,r,c);
h = pcolor(Lon,Lat,Pect);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colorbar;

subplot(2,2,2);
Pect = reshape(Pcnt_hdw_Future2,r,c);
h = pcolor(Lon,Lat,Pect);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colorbar;

subplot(2,2,3);
Pect = reshape(Pcnt_ffwi_Current2,r,c);
h = pcolor(Lon,Lat,Pect);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colorbar;

subplot(2,2,4)
Pect = reshape(Pcnt_ffwi_Future2,r,c);
h = pcolor(Lon,Lat,Pect);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colorbar;

%% CHANGE THE EXTREME DAYS

load('NumExtremeDay.mat');
load('N_HDW_CurrentCumulate.mat');
[r,c] = size(Lat);
npixel = r*c;
cumulateExtreme(:,1) = reshape(Cumulate,npixel,1);
load('N_HDW_FutureC.mat');
cumulateExtreme(:,2) = reshape(Cumulate,npixel,1);
load('N_cfwi_current.mat');
cumulateExtreme(:,3) = reshape(Cumulate,npixel,1);
load('N_cfwi_FutureC.mat');
cumulateExtreme(:,4) = reshape(Cumulate,npixel,1);
load('N_ffwi_CurrentCumulate.mat');
cumulateExtreme(:,5) = reshape(Cumulate,npixel,1);
load('N_ffwi_FutureCumulate.mat');
cumulateExtreme(:,6) = reshape(Cumulate,npixel,1);

for i = 1:6
    if (mod(i,2)==0)
        PcntBad(:,i) = cumulateExtreme(:,i)./cfwi_Length(:,2);
    else
        PcntBad(:,i) = cumulateExtreme(:,i)./cfwi_Length(:,1);
    end    
end


%% Find the 1% and 5% maximum for all indices
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Normalization\Indices\Std;
dirInfo = dir('*.mat');
load('N_cfwi_futureC.mat');
cfwi_future = TimeSeries_Stdize;
FSL_future = sum(~isnan(cfwi_future),2);
load('N_cfwi_current.mat');
cfwi_current = TimeSeries_Stdize;
FSL_current = sum(~isnan(cfwi_current),2);
load('C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2010.mat');
load('C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2100.mat');
[r,c] = size(Lat);
npixel = r*c;

LandCurrent = reshape(LandCurrent,npixel,1);
LandFuture = reshape(LandFuture,npixel,1);

for i = 1:4 % don't calculate forsberg
    load (dirInfo(i).name);
    
    % set fire season length and study region
   
    if(mod(i,2)~=0)%current HDW
        TimeSeries(isnan(cfwi_current)) = -1; 
    else % future HDW
        TimeSeries(isnan(cfwi_future)) = -1;
    end
    
    [Descend_TS,Index_D] = sort(TimeSeries,2,'descend');
    n_1percentile = fix(0.01*size(TimeSeries,2));
    percentile1_Max = Descend_TS(:,1:n_1percentile);
    mean_percentile1 = nanmean(percentile1_Max,2);

    n_5percentile = fix(0.01*size(TimeSeries,2));
    percentile5_Max = Descend_TS(:,1:n_5percentile);
    mean_percentile5 = nanmean(percentile5_Max,2);
    
    if(mod(i,2)==0)%future
         mean_percentile1(LandFuture==0)=  nan;
         mean_percentile5(LandFuture==0)=  nan;
    else
        mean_percentile1(LandCurrent==0)=  nan;
        mean_percentile5(LandCurrent==0)=  nan;
    end
    
    Pcnt1(:,i) = mean_percentile1;
    Pcnt5(:,i) = mean_percentile5;
end

DiffMax5 = Pcnt5(:,[2,4])-Pcnt5(:,[1,3]);
DiffMax1 = Pcnt1(:,[2,4])-Pcnt1(:,[1,3]);


%% Collect maximum, minimum, range
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Normalization\Indices_Separate;
load('C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\UnitedStates_Polygon.mat');
dirInfo = dir('*.mat');
n = length(dirInfo);
load ('C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2100.mat');
[r,c] = size(LandFuture);
npixel = r*c;
% figure;
T = ["HDW-Current", "HDW-Future","cfwi-Current", "cfwi-Future","ffwi-Current", "ffwi-Future"];
for i = 1 :n
    % calculate number of days > 0.5
    load(dirInfo(i).name);
    Extreme_Index = TimeSeries_Stdize>0.5;
    cumulateExtreme(:,i) = sum (Extreme_Index,2);
    Max(:,i) = reshape(Pcnt1_Max,npixel,1);
    Min(:,i) = reshape(Pcnt1_Min,npixel,1);
    Range(:,i) = reshape(range,npixel,1);
    cumulateExtreme(cumulateExtreme==0) = NaN;
end

DiffMax(:,1) = cumulateExtreme(:,2)-cumulateExtreme(:,1);
DiffMax(:,2) = cumulateExtreme(:,4)-cumulateExtreme(:,3);
DiffMax(:,3) = cumulateExtreme(:,6)-cumulateExtreme(:,5);

subplot(3,2,1);
Value = reshape(Max(:,1),r,c);
h = pcolor(Lon,Lat,Value);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
caxis([200,1200]);
title('HDW (Current)');

subplot(3,2,2);
Value = reshape(Max(:,2),r,c);
h = pcolor(Lon,Lat,Value);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
caxis([200,1200]);
title('HDW (Future)');

subplot(3,2,3);
Value = reshape(Max(:,3),r,c);
h = pcolor(Lon,Lat,Value);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
caxis([0,400]);
title('CFWI (Current)');

subplot(3,2,4);
Value = reshape(Max(:,4),r,c);
h = pcolor(Lon,Lat,Value);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
caxis([0,400]);
title('CFWI (Future)');

subplot(3,2,5);
Value = reshape(Max(:,5),r,c);
h = pcolor(Lon,Lat,Value);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
caxis([40,120]);
title('FFWI (Current)');

subplot(3,2,6);
Value = reshape(Max(:,6),r,c);
h = pcolor(Lon,Lat,Value);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
caxis([40,120]);
title('FFWI(Future)');

figure
subplot(3,1,1);
Value = reshape(DiffMax(:,1),r,c);
h = pcolor(Lon,Lat,Value);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
title('HDW (Future-Current)');

subplot(3,1,2);
Value = reshape(DiffMax(:,2),r,c);
h = pcolor(Lon,Lat,Value);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
title('CFWI (Future-Current)');

subplot(3,1,3);
Value = reshape(DiffMax(:,3),r,c);
h = pcolor(Lon,Lat,Value);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
title('FFWI (Future-Current)');

figure;
for i = 1:length(dirinfo)
    load(dirinfo(i).name);
    subplot(3,2,i);
    Value = Pect;
    h = pcolor(Lon,Lat,Value);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');
    colormap('jet');
    caxis([0,55]);
    colorbar;
end