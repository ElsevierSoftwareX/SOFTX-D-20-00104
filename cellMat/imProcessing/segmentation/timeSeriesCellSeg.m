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


function timeSeriesCellSeg(path_names,experiment_name,folder_names,cellsegParam)

% This function segments the cell for all the timepoints within the
% timelapse

f = progress_app;
pause(3);
f.UIFigure.Name = 'TFM LAB: Segmenting cells';
d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
    'Message','Creating output folders...');

%Access the files in the shift corrected folder
files = dir([path_names.store_path filesep experiment_name filesep folder_names.shift_corrected filesep folder_names.cell]);
files = files(3:end);

%Define the channel_name
channel_name = folder_names.cellSeg;

%Create a folder to store the segmentation images
save_path = [path_names.store_path filesep experiment_name filesep folder_names.shift_corrected filesep channel_name];
if not(isfolder(save_path))
    mkdir(save_path);
end

%Loop through the timepoints
for tt = 1:length(files)     
        
    %Load the timepoint
    file_name = files(tt).name;
    
    %Update progress bar
    d.Value = tt/length(files);
    d.Message = sprintf(['Segmenting cell at ' file_name(1:end-4) '...']);
    
    %Load the image
    load([files(tt).folder filesep file_name],'im_corr');
     
    % Segment the cell    
    im_corr = dip_image(im_corr);
    cell_seg = cellSegmentation(im_corr,cellsegParam);
    clear im_corr;
    cell_seg = im2mat(cell_seg>0,'bin');

    % Store the results
    save([save_path filesep file_name],'cell_seg');
    
    clear cell_seg;   

end 

% Close the dialog box
close(d);
delete(f);