function Ave2D = TemporalAvg(InputDir, OutputFileName)

    j = 0;
    for i = 1:length(InputDir)
        FileName = InputDir(i).name;
        load(FileName);    
        if (mod (str2num(FileName(end-9:end-6)),4)==0 && strcmp(FileName(end-5:end-4),'02'))

            % save data for 02/29 of each year
            if j==0
              Leap = ffwi_month(:,:,end);
            else
              Leap = cat(3,Leap,ffwi_month(:,:,end));
            end
            % clean data for 02/29 in the annually stacked dataset
            ffwi_month(:,:,end) = [];
            j = j+1;
        end

            % each file has a 259*359*364 matrix (daily stack)
        [l,w,h] = size(ffwi_month);
        % annual stack: each year is a layer
        % each column is one day for the whole region
        % each row is the annual cycle of one pixel
        month_v = reshape(ffwi_month,l*w,h);
        if (i==1)
            ffwi_month_v= month_v;
        else
            ffwi_month_v = cat(2,ffwi_month_v,month_v);
        end
        
    end
    
    AveSurface = mean(ffwi_month_v,2);
    Ave2D = reshape(AveSurface,l,w);
%     h= pcolor(Ave2D);
%     h.EdgeColor = 'none';
%     h.ZData = h.CData;
    save(OutputFileName,'Ave2D');% Seasonal average value (pixel-wise)
end

