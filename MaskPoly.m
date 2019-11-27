%%% This is a function that mask dataset based on a specific polygon %%%

function CleanData = MaskPoly (InputLat,InputLon,InputData,MaskIndex)

% % Boundary of United States
% Bndy_Lon = Polygon.Vertices(:,1);
% Bndy_Lat = Polygon.Vertices(:,2);
% 
%nIdentify whether our datasets fall within the Boundary
[n,m] = size(InputLon);
npixel = n*m;
MaskIndex_V = reshape(MaskIndex,npixel,1);
Lon_V = reshape(InputLon,npixel,1);
Lat_V = reshape(InputLat,npixel,1);
ffwi_V = reshape(InputData,npixel,1);
% % use inpolygon() to tell whether data points fall within the boundary
% MaskIndex = inpolygon(Lon_V,Lat_V,Bndy_Lon,Bndy_Lat);
% % mask out grids outside of United States: give NaN to outside grids
all = [MaskIndex_V,Lat_V,Lon_V,ffwi_V];
all((MaskIndex==0),end) = nan;
CleanData = reshape(all(:,end),n,m);
