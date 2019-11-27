cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\LandCoverType

L2010 = readtable('xyz_lnd_veg.2011_A1B.txt');
L2100 = readtable('xyz_lnd_veg.2100_A1B.txt');

Lat = L2010.lat;
Lon = L2010.lon;
LT2010 = L2010.veg;
LT2100 = L2100.veg;

Lat_Grid = reshape(Lat,359,259)';
% Lat_Grid = flipud(Lat_Grid);
Lon_Grid = reshape(Lon,359,259)';

LT2010_Grid = reshape(LT2010,359,259)';
% LT2010_Grid = flipud(LT2010_Grid);
LT2100_Grid = reshape(LT2100,359,259)';
% LT2100_Grid = flipud(LT2100_Grid);

figure

subplot(1,2,1)
h1 = pcolor(Lon_Grid,Lat_Grid,LT2010_Grid);
% hide grids
h1.EdgeColor = 'none';
% let data tips show the true value of pixels
h1.ZData = h1.CData;
title('LULC 2010');
c1 = caxis;

subplot(1,2,2)
h2 = pcolor(Lon_Grid,Lat_Grid,LT2100_Grid);
h2.EdgeColor = 'none';
h2.ZData = h2.CData;
title('LULC 2100');
c2 = caxis;

c3 = [min([c1,c2]),max([c1,c2])];
colorbar;

% colorbar('location','Manual', 'position', [0.93 0.1 0.02 0.81]);

colorbar;


% Urban 
[UrbanR,UrbanC] = find(LT2010_Grid==2);
% Cropland
[CropR,CropC] = find(LT2010_Grid==13);
% Trees
[ForestR,ForestC] = find(LT2010_Grid>7 & LT2010_Grid<11);
% Grass
[GrassR,GrassC] = find(LT2010_Grid==11);
% Shrubland
[ShrubR,ShrubC] = find(LT2010_Grid==12);

for i = 1:size(UrbanR,1)
    plot(squeeze(ffwi(UrbanR(i),UrbanC(i),:)));
    hold on;
end

for i = 1:size(ForestR,1)
    plot(squeeze(ffwi(ForestR(i),ForestC(i),:)));
    hold on;
end


imwrite(uint8(LT2010_Grid),'Current_LULC.tiff','tiff');
imwrite(uint8(LT2100_Grid),'Future_LULC.tiff','tiff');

%% write to geotiff images
[n,m] = size(Lat);
p = n*m;
Lat_V = reshape(Lat,p,1);
Lon_V = reshape(Lon,p,1);

R = georasterref('RasterSize',size(LT2010_Grid), ...
    'LatitudeLimits',[min(Lat_V) max(Lat_V)],'LongitudeLimits',[min(Lon_V) max(Lon_V)]);

 
geotiffwrite('LT2010.tif',LT2010_Grid,R);
geotiffwrite('LT2100.tif',LT2100_Grid,R);

