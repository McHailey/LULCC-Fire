%%% data format conversion %%%

%% load data %%
load('ExtremeValues.mat');
r = 259;
c = 359;
Name = ["HDW-Current","HDW-Future","CFWI-Current","CFWI-Future"];

%% to .nc file

load ('OverallResults.mat');

[y1,x1] = size(reshape(Cnsct(:,1),259,359));

% write Cnsct to nc file
nccreate('Current.nc','cfwi_Cnsct','Dimensions', {'y',y1,'x',x1});
nccreate('Current.nc','hdw_Cnsct','Dimensions', {'y',y1,'x',x1});
nccreate('Future.nc','cfwi_Cnsct','Dimensions', {'y',y1,'x',x1});
nccreate('Future.nc','hdw_Cnsct','Dimensions', {'y',y1,'x',x1});

% current cfwi
cfwi_Cnsct = reshape(Cnsct(:,1),259,359);
ncwrite('Current.nc','cfwi_Cnsct',cfwi_Cnsct);
% future cfwi
cfwi_Cnsct = reshape(Cnsct(:,2),259,359);
ncwrite('Future.nc','cfwi_Cnsct',cfwi_Cnsct);
% current hdw
hdw_Cnsct = reshape(Cnsct(:,3),259,359);
ncwrite('Current.nc','hdw_Cnsct',hdw_Cnsct);
% future hdw
hdw_Cnsct = reshape(Cnsct(:,4),259,359);
ncwrite('Future.nc','hdw_Cnsct',hdw_Cnsct);

% write Percentage of Bad Days to nc files
nccreate('Current.nc','cfwi_PcntBad','Dimensions', {'y',y1,'x',x1});
nccreate('Current.nc','hdw_PcntBad','Dimensions', {'y',y1,'x',x1});
nccreate('Future.nc','cfwi_PcntBad','Dimensions', {'y',y1,'x',x1});
nccreate('Future.nc','hdw_PcntBad','Dimensions', {'y',y1,'x',x1});

% current cfwi
cfwi_PcntBad = reshape(PcntBad(:,1),259,359);
ncwrite('Current.nc','cfwi_PcntBad',cfwi_PcntBad);
% future cfwi
cfwi_PcntBad = reshape(PcntBad(:,2),259,359);
ncwrite('Future.nc','cfwi_PcntBad',cfwi_PcntBad);
% current hdw
hdw_PcntBad = reshape(PcntBad(:,3),259,359);
ncwrite('Current.nc','hdw_PcntBad',hdw_PcntBad);
% future hdw
hdw_PcntBad = reshape(PcntBad(:,4),259,359);
ncwrite('Future.nc','hdw_PcntBad',hdw_PcntBad);

% write Averaged 5% Max to nc file
nccreate('Current.nc','cfwi_MaxPcnt5','Dimensions', {'y',y1,'x',x1});
nccreate('Current.nc','hdw_MaxPcnt5','Dimensions', {'y',y1,'x',x1});
nccreate('Future.nc','cfwi_MaxPcnt5','Dimensions', {'y',y1,'x',x1});
nccreate('Future.nc','hdw_MaxPcnt5','Dimensions', {'y',y1,'x',x1});

% current cfwi
cfwi_MaxPcnt5 = reshape(AvgMaxPcnt5(:,1),259,359);
ncwrite('Current.nc','cfwi_MaxPcnt5',cfwi_MaxPcnt5);
% future cfwi
cfwi_MaxPcnt5 = reshape(AvgMaxPcnt5(:,2),259,359);
ncwrite('Future.nc','cfwi_MaxPcnt5',cfwi_MaxPcnt5);
% current hdw
hdw_MaxPcnt5 = reshape(AvgMaxPcnt5(:,3),259,359);
ncwrite('Current.nc','hdw_MaxPcnt5',hdw_MaxPcnt5);
% future hdw
hdw_MaxPcnt5 = reshape(AvgMaxPcnt5(:,4),259,359);
ncwrite('Future.nc','hdw_MaxPcnt5',hdw_MaxPcnt5);

% write p-value
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices\Analysis
load ('ttest.mat');
nccreate('PValue.nc','p_cfwi','Dimensions', {'y',y1,'x',x1});
nccreate('PValue.nc','p_hdw','Dimensions', {'y',y1,'x',x1});

p_cfwi = reshape(p_cfwi,259,359);
ncwrite('PValue.nc','p_cfwi',p_cfwi);

p_hdw = reshape(p_hdw,259,359);
ncwrite('PValue.nc','p_hdw',p_hdw);



%% to csv %%
% for i = 1:size(AvgExtreme,2)
%     output = reshape(AvgExtreme(:,i),r,c);
%     writematrix(output,strcat(Name(i),".csv"),'Delimiter','tab')
% end
% Output = [Lat_V,Lon_V,AvgExtreme];
% writematrix(Output,"NumBadDays.csv",'Delimiter','tab');
% writematrix(Lat,"Latitude.csv",'Delimiter','tab');
% writematrix(Lon,"Longitude.csv",'Delimiter','tab');

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\LandCoverType

L2010 = readtable('xyz_lnd_veg.2011_A1B.txt');
L2100 = readtable('xyz_lnd_veg.2100_A1B.txt');

r = 259;
c = 359;
npixel = r*c;

% add variables to LULCC Table (L2010/L2100)
L2010.cfwi_Cnsct = reshape(reshape(Cnsct(:,1),r,c)',npixel,1);
L2010.HDW_Cnsct = reshape(reshape(Cnsct(:,3),r,c)',npixel,1);
L2010.cfwi_PcntBad = reshape(reshape(PcntBad(:,1),r,c)',npixel,1);
L2010.HDW_PcntBad = reshape(reshape(PcntBad(:,3),r,c)',npixel,1);
L2010.cfwi_Max5 = reshape(reshape(AvgMaxPcnt5(:,1),r,c)',npixel,1);
L2010.HDW_Max5 = reshape(reshape(AvgMaxPcnt5(:,3),r,c)',npixel,1);


L2100.cfwi_Cnsct = reshape(reshape(Cnsct(:,2),r,c)',npixel,1);
L2100.HDW_Cnsct = reshape(reshape(Cnsct(:,4),r,c)',npixel,1);
L2100.cfwi_PcntBad = reshape(reshape(PcntBad(:,2),r,c)',npixel,1);
L2100.HDW_PcntBad = reshape(reshape(PcntBad(:,4),r,c)',npixel,1);
L2100.cfwi_Max5 = reshape(reshape(AvgMaxPcnt5(:,2),r,c)',npixel,1);
L2100.HDW_Max5 = reshape(reshape(AvgMaxPcnt5(:,4),r,c)',npixel,1);

writetable (L2010,'CurrentFull.txt');
writetable (L2100,'FutureFull.txt');

% save another file for ttest
r = 259;
c = 359;
npixel = r*c;
i = L2010.i;
j = L2010.j;
lat = L2010.lat;
lon = L2010.lon;
p_cfwi = reshape(reshape(p_cfwi,r,c)',npixel,1);
p_hdw = reshape(reshape(p_hdw,r,c)',npixel,1);
t = table(i,j,lat,lon,p_cfwi,p_hdw);

writetable (t,'ttestMax5%.txt');

%% to tif %%
p = r*c;
Lat_V = reshape(Lat,p,1);
Lon_V = reshape(Lon,p,1);

R = georasterref('RasterSize',size(Lat), ...
    'LatitudeLimits',[min(Lat_V) max(Lat_V)],'LongitudeLimits',[min(Lon_V) max(Lon_V)]);

for i = 2:size(AvgExtreme,2)
    output = reshape(AvgExtreme(:,i),r,c);
    geotiffwrite(strcat(Name(i),".tif"),output,R);
end

geotiffwrite("Latitude.tif",Lat,R);
geotiffwrite("Longitude.tif",Lon,R);


% h = pcolor(Lon,Lat,reshape(AvgExtreme(:,1),r,c));
% h.EdgeColor = 'none';
% h.ZData =h.CData; 
% hold on;
% plot([states.Lon],[states.Lat], 'Color','black');
% colormap('jet');
% colorbar;
% title('ffwi (future), Frequency of Large Fire Days within CFWI Fire Season');

%% To plots %%

% Percentage of Bad Days

cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices\Analysis
load('ExtremeValue.mat');
PcntBad(:,[1,3]) = 100*Bad(:,[1,3])./cfwi_FL(:,1);%current
PcntBad(:,[2,4]) = 100*Bad(:,[2,4])./cfwi_FL(:,2);%future

DiffPcntBad = PcntBad(:,[2,4])-PcntBad(:,[1,3]);

% DiffMax5 = Pcnt5(:,[2,4])-Pcnt5(:,[1,3]);
% DiffNumBad = AvgNumBad(:,[2,4,6])-AvgNumBad(:,[1,3,5]);
% DiffConsecutive = AvgConsecutive(:,[2,4])-AvgConsecutive(:,[1,3]);
% ExtreDiff = ExtremeNum(:,[2,4])-ExtremeNum(:,[1,3]);

j = 0;
figure;
for i = 1:6
  if (mod(i,3)==0)
    subplot(2,3,i);
    Value = reshape(DiffPcntBad(:,i/3),r,c);
    h = pcolor(Lon,Lat,Value);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');
    
%     % caxis for maximum
%     if(i==3)
%         caxis ([-200,200]);
%     elseif(i==6)
%         caxis([-50,50]);
%     end
%     caxis([-100,100]);
    caxis([-10,10]);
    colormap(gca, bluewhitered);
    colorbar;
    xlim([-130,-65]);
    ylim([25,50]);
    title("Future-Current");
      
  else
    j = j+1;  
    subplot(2,3,i);
    Value = reshape(PcntBad(:,j),r,c);
    h = pcolor(Lon,Lat,Value);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');
   
%     caxis([50,500]);
    xlim([-130,-65]);
    ylim([25,50]);
%       % caxis for maximum
%     if (j <=2)
%         caxis([0,1000])
%     elseif(j>2&&j<=4)
%         caxis([0,300]);
%     end
%     if(j<=2)
%         caxis([0,100]);
%     else
%         caxis([0,120]);
%     end
    caxis([0,55]);
    colormap(gca, 'jet');
    colorbar;
    title(T(j));
  end
end

% consecutive bad days
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices\Analysis
dirCnsct = dir('Cnsct*.mat');


for i = 1:length(dirCnsct)
    load(dirCnsct(i).name);
    Cnsct(:,i) = nCnsct;
    T(i) = string(dirCnsct(i).name(9:end-4));
end

AvgCnsct = Cnsct./10;
CnsctDiff = AvgCnsct(:,[2,4])-AvgCnsct(:,[1,3]);


j=0;
figure;
for i = 1:6
  if (mod(i,3)==0)
    subplot(2,3,i);
    Value = reshape(CnsctDiff(:,i/3),r,c);
    h = pcolor(Lon,Lat,Value);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');
    
%     % caxis for maximum
%     if(i==3)
%         caxis ([-200,200]);
%     elseif(i==6)
%         caxis([-50,50]);
%     end
%     caxis([-100,100]);
    caxis([-5,5]);
    colormap(gca, bluewhitered);
    colorbar;
    xlim([-130,-65]);
    ylim([25,50]);
    title("Future-Current");
      
  else
    j = j+1;  
    subplot(2,3,i);
    Value = reshape(AvgCnsct(:,j),r,c);
    h = pcolor(Lon,Lat,Value);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');

    xlim([-130,-65]);
    ylim([25,50]);
% for average number of consecutive bad days
    caxis([0,15]);
    colormap(gca, 'jet');
    colorbar;
    title(T(j));
  end
end

% 5% maximum averaged value
cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireIndices\Analysis
load('ExtremeValue.mat');

Max5Diff = 100.*(AvgMaxPcnt5(:,[2,4])-AvgMaxPcnt5(:,[1,3]))./AvgMaxPcnt5(:,[1,3]);

j=0;
figure;
for i = 1:6
  if (mod(i,3)==0)
    subplot(2,3,i);
    Value = reshape(Max5Diff(:,i/3),r,c);
    h = pcolor(Lon,Lat,Value);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');
    
    % caxis for relative changes of 5% maximum value
    % (Future-current)/current
    caxis ([-30,30]);

    colormap(gca, bluewhitered);
    colorbar;
    xlim([-130,-65]);
    ylim([25,50]);
    title("100*(Future-Current)/Current");
      
  else
    j = j+1;  
    subplot(2,3,i);
    Value = reshape(AvgMaxPcnt5(:,j),r,c);
    h = pcolor(Lon,Lat,Value);
    h.EdgeColor = 'none';
    % h4.ZData =h4.CData; 
    hold on;
    plot([states.Lon],[states.Lat], 'Color','black');

    xlim([-130,-65]);
    ylim([25,50]);
    
    % caxis for maximum
    if(j<=2)
        caxis([0,250]);
    else
        caxis([0,1000]);
    end

    colormap(gca, 'jet');
    colorbar;
    title(T(j));
  end
end