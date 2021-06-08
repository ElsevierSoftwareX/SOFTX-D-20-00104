%    This file is part of TFMLAB.
%     Copyright (C) 2016-2020 Bme, Dep. Mech. Engineering, KUleuven (Belgium)
%     Copyright (C) 2016 Alvaro Jorge-Penas
%     Copyright (C) 2019-2020 Jorge Barrasa-Fano
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


function timeSeriesImFilt(store_path,experiment_name,folder_names,file_names,filtParam,cellsegParam,cellfiltParam)
% This function applies a filtering step to all the timepoints within the
% timelapse
f = progress_app;
pause(3);
f.UIFigure.Name = 'TFM LAB: Filtering';
d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
    'Message','Creating output folders...');
%Create the filtering results folder and store the paramtere file
save_path = [store_path filesep experiment_name filesep folder_names.tiff_filtered];
if not(isfolder(save_path))
    mkdir(save_path);
end
save([save_path filesep file_names.filt_param],'filtParam','cellsegParam','cellfiltParam');

%Check which channels we have
folders = dir([store_path filesep experiment_name filesep folder_names.tiff_raw]);
folders = folders(3:end);

%Loop through the channels
for ch = 1:length(folders)
    
    if folders(ch).isdir
        
        %Get the channel name
        channel_name = folders(ch).name;
        
        %Get the files within this folder
        files = dir([folders(ch).folder filesep folders(ch).name]);
        files = files(3:end);
        
        %Define an auxiliary counter for the progress bar
        cnt = 1;
        
        %Loop through the time points
        for tt = 1:length(files)
            %Get the file name
            file_name = files(tt).name;
            
            %Update progress bar
            d.Value = cnt/(length(folders)*length(files));
            d.Message = sprintf(['Filtering ' channel_name ' at ' file_name(1:end-4) '...']);
            cnt = cnt +1;
            
            %Load the image
            im = read3Dim([files(tt).folder filesep file_name]);
            switch channel_name
                case folder_names.beads
                    im = beadFilter(im,filtParam.beads.var1,filtParam.beads.var2,filtParam.beads.weight,filtParam.beads.smoothVar,filtParam.beads.stretchLB,filtParam.beads.stretchUB);
                case folder_names.cell
                    [im,cellfiltParam.smoothVal] = cellFilter(im,cellfiltParam.lightenDarkRegions,cellfiltParam.smoothVal,cellfiltParam.stretchLB,cellfiltParam.stretchUB);
                case folder_names.fibers
                    [im,~] = fiberFilter(im,filtParam.fibers.lightenDarkRegions,filtParam.fibers.smoothVal,filtParam.fibers.stretchLB,filtParam.fibers.stretchUB);
            end
            %Create a saving folder
            tmp_save_path = [save_path filesep channel_name];
            if not(isfolder(tmp_save_path))
                mkdir(tmp_save_path);
            end
            save([tmp_save_path filesep file_name(1:end-4) '.mat'],'im');
            clear tmp_save_path im;
        end
    end
end
% Close the dialog box
close(d);
delete(f);