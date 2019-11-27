function ReadData(FolderPath,FireIndex)

% cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2010\Raw

cd(FolderPath);
dinfo = dir('*.nc');

for i = 1:size(dinfo,1)
    Index(:,:,i) = ncread(dinfo(i,:).name,FireIndex)';
    
%   HDW doesn't have Lat/Lon
%   Lat = ncread(dinfo(1,:).name,'XLAT')';
%   Lon = ncread(dinfo(1,:).name,'XLONG')';
%   h = pcolor(Lon,Lat,Index(:,:,i));
%   h.EdgeColor = 'none';
%   h.ZData = h.CData;
end

StartYear = 2006; %start year of dataset
EndYear = 2015;%end year of dataset
Month1 = 1;
Month2 = 12;
Day1 =1;

i = 0;
N = nan(10,1);
DayNum = nan(10,1);
NDay = 0;

for Year = StartYear:EndYear
  mDay = 0;
  i = i+1;
  if mod(Year,4)==0
      nDay = 361;
  else
      nDay = 360;
  end
  DayNum(i,1)=nDay;
  NDay = NDay+nDay;
  disp(NDay);
  Index_year = Index(:,:,NDay-nDay+1:NDay);
  % Group Data by year
  save([FireIndex,num2str(Year)],'Index_year');
%   load(['ffwi' num2str(Year) '.mat']);

% % Day1/Day2 for cfwi
%  if strcmp(FireIndex,'cfwi')
%       for Month = Month1: Month2
%       
%         if(Month==1)
%             Day2 = 28;
%         elseif (Month==2)
%             if mod(Year,4) == 0 %in the first year, start in the start month rather than January
%                 Day2 = 29;
%             else
%                 Day2=28;
%             end
%         elseif (Month==12)
%             Day2 = 29;
%         elseif ((mod(Month,2)==0 && Month>7) ||(mod(Month,2)==1 && Month<=7))
%             Day2 = 31;
%         else
%             Day2 = 30;
%         end
%         
%         mDay = mDay + Day2;
%         Index_month = Index_year(:,:,(mDay-Day2+1):mDay);
%         save([FireIndex,num2str(Year), num2str(Month,'%02d')],'Index_month');
%       end
%  % Day1/Day2 for other fire indices   
%  else
%     for Month = Month1: Month2
%         if (Month==2)
%             if mod(Year,4) == 0 %in the first year, start in the start month rather than January
%                 Day2 = 29;
%             else
%                 Day2=28;
%             end
%         % in this case, the model output doesn't produce 12/31 data 
%         elseif ((mod(Month,2)==0 && Month>7 && Month<12 ) ||(mod(Month,2)==1 && Month<=7))
%             Day2 = 31;
%         else
%             Day2 = 30;
%         end
%         
%         mDay = mDay + Day2;
%         Index_month = Index_year(:,:,(mDay-Day2+1):mDay);
%         save([FireIndex,num2str(Year), num2str(Month,'%02d')],'Index_month');
%     end
%  end
end


