%%% normalize future dataset with current ranges %%%

load ('N_ffwi_Future.mat');
FutureTS = TimeSeries;
FutureNTS = TimeSeries_Stdize;

load ('N_ffwi_Current.mat');
[r,c] = size(Lat);
npixel = c*r;
range_V = reshape(range, npixel,1);
Pcnt1_MinV = reshape(Pcnt1_Min,npixel,1);
TimeSeries = FutureTS;
TimeSeries_Stdize = (FutureTS-Pcnt1_MinV)./range_V;
Avg10 = reshape (nanmean(TimeSeries_Stdize,2),259,359);

save('N_ffwi_futureC.mat','Avg10','TimeSeries_stdize','Lat','Lon','TimeSeries');

load ('C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\UnitedStates_Polygon.mat');


load ('N_cfwi_current.mat');
cfwi = TimeSeries_Stdize;
% load('N_ffwi_FutureC.mat');
load('N_ffwi_current.mat');
HDW = 100.*TimeSeries_Stdize;
HDW(isnan(cfwi)) = nan;
Extreme_Index= HDW>50;
% number of large events (>50%), future normalized by current
Cumulate_Index= nansum(Extreme_Index,2);
Cumulate_Index(Cumulate_Index==0)= nan;

load('NumExtremeDay.mat');
% percentage of large events based on CFWI fire season length
Pcnt_cfwi_Future2 = 100.*Cumulate_Index./cfwi_Length(:,2);

Cumulate2 = reshape(Cumulate_Index,r,c);
Pect = reshape(Pcnt_cfwi_Future2,r,c);

h = pcolor(Lon,Lat,Cumulate2);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
title('ffwi (future), Frequency of Large Fire Days within CFWI Fire Season');

