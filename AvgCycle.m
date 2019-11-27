%%% use dataset grouped by year to calculate the mean annual cycle

function AvgOutput = AvgCycle(InputFolder,OutputFileName)

% for 2/29, only leap years have it
x = [2019 2 29];
x_date = datetime(x);
x_doy = day(x_date,'dayofyear');

cd(InputFolder);
% cd C:\Users\meicheng\Desktop\PhD-year1\PhysicalGeography\Project\FireData2010\Annual
dirInfo = dir('*.mat');
j = 0;
for i = 1:length(dirInfo)
    FileName = dirInfo(i).name;
    load(FileName);    
    if (mod (str2num(FileName(end-7:end-4)),4)==0)
        
        % save data for 02/29 of each year
        if j==0
          Leap = ffwi_year(:,:,x_doy);
        else
          Leap = cat(3,Leap,ffwi_year(:,:,x_doy));
        end
        % clean data for 02/29 in the annually stacked dataset
        ffwi_year(:,:,x_doy) = [];
        j = j+1;
    end
    
        % each file has a 259*359*364 matrix (daily stack)
    [l,w,h] = size(ffwi_year);
    % annual stack: each year is a layer
    % each column is one day for the whole region
    % each row is the annual cycle of one pixel
    ffwi_year_v(:,:,i) = reshape(ffwi_year,l*w,h);
end

% calculate the average daily annual cycle (without 02/29)
AvgOutput = mean(ffwi_year_v,3);
save(OutputFileName,'AvgAnnualCycleOutput');

end