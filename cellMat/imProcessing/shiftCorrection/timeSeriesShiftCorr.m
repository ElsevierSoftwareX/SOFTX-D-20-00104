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

function timeSeriesShiftCorr(path_names,experiment_name,folder_names,shift_info,rot_info)
% This function applies a given shift (rigid translation) and a rotation to each timepoint of the timelapse

f = progress_app;
pause(3);
f.UIFigure.Name = 'TFM LAB: Shift correction';
d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
    'Message','Creating output folders...');

%Access the files in the folder
folders = dir([path_names.store_path filesep experiment_name filesep folder_names.tiff_filtered]);
folders = folders(3:end);

%Keep only the ones that are actually folders 
indx_rmv = [];
for ii = 1:length(folders)
    if not(folders(ii).isdir)
        indx_rmv = [indx_rmv ii];
    end
end
folders(indx_rmv) = [];
%Compute the number of time points
n_channels = length(folders);

% Compute the total required crop
cropStart = zeros(size(shift_info.ref));
cropEnd = zeros(size(shift_info.ref));

cropStart(shift_info.ref>0) = ceil(shift_info.ref(shift_info.ref>0));
cropEnd(shift_info.ref<0) = ceil(abs(shift_info.ref(shift_info.ref<0)));

cropStartT = max(cropStart,[],1) ;
cropEndT = max(cropEnd,[],1);

%Create a folder to store the shift corrected images and store the
%parameter file
save_path = [path_names.store_path filesep experiment_name filesep folder_names.shift_corrected];
if not(isfolder(save_path))
    mkdir(save_path);
end


for ch = 1:n_channels    
    
    %Get the channel name
    channel_name = folders(ch).name;
    
    %Get the files within the folder
    files = dir([folders(ch).folder filesep channel_name]);
    files = files(3:end);
    nTP = length(files);
    
    %Loop through the time points to get 
    for tt = 1:nTP
        
        %Load the timepoint
        file_name = files(tt).name;
        load([files(tt).folder filesep file_name],'im');      
        
        %Update progress bar
        d.Value = tt/nTP;
        d.Message = sprintf(['Correcting shifts for channel ' channel_name ' ' file_name(1:end-4) '...']);
        
        % Correct the shift     
        im_corr = shiftCorr(im,shift_info.ref(tt,:),cropStartT,cropEndT,'im'); 
        
        %Clear the image to save RAM
        clear im;
        
        %Create a saving folder
        tmp_save_path = [save_path filesep channel_name];
        if not(isfolder(tmp_save_path))
            mkdir(tmp_save_path);
        end
        save([tmp_save_path filesep file_name],'im_corr');
        
        clear im_corr;  

    end
    
end

% Close the dialog box
close(d);
delete(f);



