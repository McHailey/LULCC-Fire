%% find consecutive fire days (>5 Days) %%
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Normalization\Indices;

% load('N_hdw_current.mat');
load('N_HDW_futureC.mat');
ffwi = TimeSeries_Stdize;
load('N_cfwi_futureC.mat');
% load('N_cfwi_current.mat');
cfwi = TimeSeries_Stdize;

ffwi(isnan(cfwi)==1) = nan;
Extreme_Index = ffwi>0.5;

[r,c] = size(Lat);
npixel = r*c;

load('C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\Mask2100.mat');
LandFuture = reshape(LandFuture,npixel,1);

% diff([0 find(diff(Extreme_Index)) numel(Extreme_Index)])

% only works for row vectors
% d = [true, diff(Extreme_Index) ~= 0, true];  % TRUE if values change
% n = diff(find(d));               % Number of repetitions
% Y = repelem(n, n);

[RIndex,CIndex]= find(Extreme_Index==1);
Index = sortrows([RIndex,CIndex],1);% RowIndex-pixel. ColumnIndex-TimeSeries/Day
Diff = diff(Index);
Index = [Index,[nan(1,2);Diff]];

pixel = unique(Index(:,1));
ncount = nan(length(pixel),3600);
for i = 1:length(pixel)

    if mod(i,1000)==0
        disp(i);
    end
    
    if(isnan(LandFuture(i)))
        ncount(i,:) = nan;
    else
        subIndex = Index(Index(:,1)==pixel(i),:);
        nday = 1;
        ncount(i,1) = pixel(i);% location/pixel(row Index)
        for j = 1:size(subIndex,1)
            if (subIndex(j,3)==0 && subIndex(j,4)==1)
                nday = nday+1;
            else
                nday = 1;
            end
            ncount(i,j+1) = nday;
        end 
    end
end

nEvent = sum((ncount==5),2);
nEvent = [ncount(:,1),nEvent];

[r,c] = size(Lat);
npixel = r*c;
FullPixel = nan(npixel,2);
FullPixel(:,1) = (1:npixel)';
for i = 1:length(nEvent)
    FullPixel(FullPixel(:,1)==nEvent(i,1),2) = nEvent(i,2);
end

nLasting = reshape(FullPixel(:,2),r,c);

h = pcolor(Lon,Lat,nLasting);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
colormap('jet');
colorbar;
title('Number of Lasting Large Events among 10 Years (ffwi,future)');

% differnces
load('LastingEvent_FutureC.mat');
future = nLasting;
load('LastingEvent_Current.mat');
current = nLasting;
deltaFreq2 = future - current;
h = pcolor(Lon,Lat,deltaFreq2);
h.EdgeColor = 'none';
% h4.ZData =h4.CData; 
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
colormap(bluewhitered);
colorbar;
title('Differences of Lasting Large Fire Events (Future-Current)');


%% put all consecutive days into one matrix

load('hdw_Lasting_Current.mat');
[r,c] = size(nLasting);
npixel = r*c;
consecutive(:,1) = reshape(nLasting,npixel,1);
load('hdw_Lasting_Future.mat');
consecutive(:,2) = reshape(nLasting,npixel,1);
load('ffwi_Lasting_Current.mat');
consecutive(:,3) = reshape(nLasting,npixel,1);
load('ffwi_Lasting_Future.mat');
consecutive(:,4) = reshape(nLasting,npixel,1);

AvgConsecutive = ceil(consecutive./10);

