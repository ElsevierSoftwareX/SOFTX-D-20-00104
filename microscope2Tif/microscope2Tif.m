%    This file is part of TFMLAB.
%     Copyright (C) 2020 MAtrix
%     Copyright (C) 2020 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2020 Esc. Téc. Sup. de Ingeniería, Universidad de Sevilla (Spain)
%
%     This library is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published
%     by the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This software is provided "as is",
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details
%     <http://www.gnu.org/licenses/>.




function folder_names = microscope2Tif(file_path,input_file_name,store_path,folder_names,param,file_names)
     
f = progress_app;
pause(3);
f.UIFigure.Name = 'TFM LAB: Input definition';
d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
    'Message','Creating output folders...');


% Open the bioformats reader from the input file
reader = bfGetReader([file_path filesep input_file_name]);

%Store the number of z planes
zPlanes = param.zPlanes;
tPoints = param.tPoints;
%Check if there is a relaxed separated series
if param.series_id_relaxed >= 0 
    %There is a separated relaxed file
    series_ids = [param.series_id_stressed param.series_id_relaxed];
else
    %The relaxed file is within the stressed file
    series_ids = param.series_id_stressed;
end

%Define an auxiliary container for the number of time points in each series
tPoints_container = zeros(1,length(series_ids));

% Loop on Series to get metadata
for ss = 1:length(series_ids) 
    
    % Select the current series
    info = extractMetadata(reader,series_ids(ss));
    
    %Get the Z planes to be used
    zPlanes = param.zPlanes(param.zPlanes<=info.nZplanes);
    info.nZplanes = length(zPlanes);    
    
    %Get the timepoints to be used
    tPoints = param.tPoints(param.tPoints<=info.nTsteps);
    
    %Get the x and y ranges to be used
    xRange = param.xRange(param.xRange<=info.imSize(2));        
    yRange = param.yRange(param.yRange<=info.imSize(1));
        
    %Save the number of time points that this series has:
    tPoints_container(ss) = tPoints(end);
    
    % Get the image resolution
    metadataStr = getResMetadata(info);   
    
    % Loop on time points
    for tt = tPoints
        
%         disp(['Exporting time point ' num2str(tt)]);
        
        %Update progress bar
        d.Value = tt/length(tPoints);
        d.Message = sprintf(['Series ' num2str(ss) ': Exporting time point ' num2str(tt) '...']);

        % Loop on channels
        for cc=1:length(param.channels_id)
            
            %Get the channel name
            channel_name = param.channels_name{cc};
            
            d.Message = sprintf(['Series ' num2str(ss) ': Exporting time point ' num2str(tt) '(' channel_name ' channel)...']);
                        
            %Define the output folder
            save_path = [store_path filesep param.experiment_name filesep folder_names.tiff_raw metadataStr filesep channel_name];
            if not(isfolder(save_path))
                mkdir(save_path);
            end
            
            %Check if we are in the stressed or in the relaxed state
            if ( (length(series_ids) == 2) && (ss == 1) ) ||...
                    ( (length(series_ids) == 1) && tt<tPoints(end) )
                %This series is the stressed state
                if tt<10
                    output_file_name = [save_path filesep 'tp_0' num2str(tt) '.' 'tif'];
                else
                    output_file_name = [save_path filesep 'tp_' num2str(tt) '.' 'tif'];
                end
            else
                %This series is the relaxed state
                if tt<length(tPoints)
                    %This relaxed state has multiple time points. We will use them as stressed points by continuing the previous count.
                    %We need to check the previous series' number of time
                    %points (tPoints_container) to continue counting.
                    if (tPoints_container(ss-1) + tt)<10
                        output_file_name = [save_path filesep 'tp_0' num2str(tPoints_container(ss-1) + tt) '.' 'tif'];
                    else
                        output_file_name = [save_path filesep 'tp_' num2str(tPoints_container(ss-1) + tt) '.' 'tif'];
                    end
                else
                    
                    output_file_name = [save_path filesep 'tp_R.' 'tif'];
                    
                end
            end
            
            % Loop on Zslices
            for zz=zPlanes % Loop on the Zstack
                
                % Read the current image
                im = readHyperstack(reader,info,param.channels_id(cc),zz,tt);
                
                %Keep only the desired range
                im = im(yRange,xRange);
                
                % Write the output file with the modified name
                if zz==param.zPlanes(1)
                    imwrite(im,output_file_name);
                else
                    imwrite(im,output_file_name,'WriteMode','append');
                end
                clear im;
                
            end % end of the zz loop
            
            d.Message = sprintf(['Series ' num2str(ss) ': Finished ' channel_name ' channel...']);
%             disp(['... Finished ' channel_name ' channel...']);
        end % end of the cc loop        
        d.Message = sprintf(['Series ' num2str(ss) ': Finished time point ' num2str(tt) '!']);
%         disp(['... Finished time point ' num2str(tt) '!']);
    end % end of the tt loop
    
    clear info;
    
end % end of the ss loop

% close the reader
reader.close();
clear reader;

%Return the new folder name
folder_names.tiff_raw = [folder_names.tiff_raw metadataStr];

%Save the param file
save([store_path filesep param.experiment_name filesep folder_names.tiff_raw filesep file_names.input_param],'param');


% Close the dialog box
close(d);
delete(f);

end % end of the main function