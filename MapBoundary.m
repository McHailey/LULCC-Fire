%%% This script is used to show the USA map %%%
%% -------------The United States Map within MATLAB package-------------------
figure
ax = usamap('conus');
states = shaperead('usastatelo', 'UseGeoCoords', true,...
  'Selector',...
  {@(name) ~any(strcmp(name,{'Alaska','Hawaii'})), 'Name'});
geoshow( states, 'DisplayType', 'polygon');
hold on;
h = pcolor(Lon,Lat,ffwi);
h.EdgeColor = 'none';

% try to union all states 
for i = 1:length(states)
    poly(i,1) = polyshape(states(i).Lon,states(i).Lat);
end
polyout = union(poly);

% the whole united states
figure
plot(polyout);

%% ----------------Clipping (wrapped into function--MapPolygon)------------------
% Boundary of United States
Bndy_Lon = polyout.Vertices(:,1);
Bndy_Lat = polyout.Vertices(:,2);

% Identify whether our datasets fall within the Boundary
[n,m] = size(Lon);
npixel = n*m;
Lon_V = reshape(Lon,npixel,1);
Lat_V = reshape(Lat,npixel,1);
ffwi_V = reshape(ffwi,npixel,1);
% use inpolygon() to tell whether data points fall within the boundary
in = inpolygon(Lon_V,Lat_V,Bndy_Lon,Bndy_Lat);
% mask out grids outside of United States: give NaN to outside grids
all = [in,Lat_V,Lon_V,ffwi_V];
all((in==0),end) = nan;
ffwi_US = reshape(all(:,end),n,m);
MaskIndex =reshape(in,n,m);

% --------------------------------------------------------------------------
% % load US boundary
% load('UnitedStates_Polygon.mat');
% load MaskIndex

% load Our fire index dataset with latitude and longitude

ffwi_US = MaskPoly (Lat,Lon,ffwi,MaskIndex);

% ----------------------------mapping---------------------------------
contourf(Lon,Lat,ffwi_US,'LineStyle','none');
% h = pcolor(Lon,Lat,ffwi_US);
% h.ZData = h.CData;
% h.EdgeColor = 'none';
hold on;
plot(polyout);
title("ffwi, US, 20060101");
hold on;
plot([states.Lon],[states.Lat], 'Color','black');
xlim([-130,-65]);
ylim([23,52]);

figure;
plot(polyout);

figure;
h = pcolor(Lon,Lat,ffwi);
h.EdgeColor = 'none';
