


%Load the displacement field
load([files(tt).folder filesep files(tt).name],'dispField');

%Load the crop_rect, the pad_size and the sampling_indx


%Define the channel_name
channel_name =  'beads';

%Define the resolution
res = [0.57 0.57 0.5];

%Define the name of the file where the results will be stored (ends with
%.vtk): use the full path
save_file = ['D:\cellMat2019_testing\test_janne3\disp_' channel_name '_tp_01.vtk'];

%...




%Crop the displacement field
dispField.X = dispField.X(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
dispField.Y = dispField.Y(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
if flag_3D
    dispField.Z = dispField.Z(crop_rect(2):crop_rect(2)+crop_rect(4)-1,crop_rect(1):crop_rect(1)+crop_rect(3)-1,:);
end

%Add the pad
dispField.X = padarray(dispField.X,pad_size*ones(1,length(res)),'replicate');
dispField.Y = padarray(dispField.Y,pad_size*ones(1,length(res)),'replicate');
if flag_3D
    dispField.Z = padarray(dispField.Z,pad_size*ones(1,length(res)),'replicate');
end

disp2vtk(dispField,res,sampling_indx,channel_name,save_file)