%Run the deafult parameters
parameterFile;

%% CHECK IF SUITESPARSE HAS BEEN COMPILED
if not(isfile([pwd filesep code_names.external filesep 'SuiteSparse' filesep 'SPQR' filesep 'MATLAB' filesep 'spqr_solve.mexw64']))
    
    %Check if there is a C++ compiler
    addons = matlab.addons.installedAddons;
    if (sum(ismember(addons.Identifier,'ML_MINGW')))        
        disp('Installing TFMLAB: SUITESPARSE...')
        current_path = pwd;
        cd external/SuiteSparse/SPQR/MATLAB
        spqr_install
        cd(current_path)
    else 
        disp('Please install the MATLAB Add-On MinGW: Add-Ons?Get Add-Ons and search for MinGW ');
    end
else
    disp('There is nothing to install' );
end

disp('TFMLAB is ready to use');

%% Check if DIPImage exists in the computer
% 
% %Check if there is a file where the path of DIPImage is stored
% if not(isfile([pwd filesep code_names.utilities filesep file_names.dip_image_path_file]))
%     %The file is empty, we need to look for it.
%     
%     f = progress_app;
%     f.UIFigure.Name = 'TFMLAB: checking libraries';
%     d = uiprogressdlg(f.UIFigure,'Title','Please Wait',...
%         'Message','Looking for the DIPImage library...','Indeterminate','on');
%     %Look for it
%     path_file = checkDipImage();
%     % Close the dialog box
%     close(d);
%     delete(f);
%     
%     %If we don't find it, we hav to ask the user where it is or if she
%     %wants to download it
%     if isempty(path_file)
%         %Asking the user
%         f = progress_app;
%         f.UIFigure.Name = 'TFMLAB: checking libraries';
%         
%         opts = {'Download','Select start file location','Cancel'};
%         selection = uiconfirm(f.UIFigure,compose('The DIPimage library was not found. You can: \n - Specify the location of the DIPimage start file. \n - Download the latest version, install and run TFMLAB again.'),'DIPImage not found','Icon','Warning','Options',opts);
%         
%         switch selection
%             case opts{1}
%                 web('http://www.diplib.org/download','-browser');
%                 close(f.UIFigure);
%                 return;
%             case opts{2}
%                 [file,path] = uigetfile('../*.m','Please select the file dipstart.m in your DIPimage folder');
%                 path_file = [path file];
%                 close(f.UIFigure);
%                 %We found dip_image, we add it to the path by running the startup file of dip_image
%                 run(path_file);
%                 %We save the path for later runs
%                 save([pwd filesep code_names.utilities filesep file_names.dip_image_path_file],'path_file');
%             case opts{3}
%                 close(f.UIFigure);
%                 return;
%         end
%     else
%         %We found dip_image, we add it to the path by running the startup file of dip_image
%         run(path_file);
%         %We save the path for later runs
%         save([pwd filesep code_names.utilities filesep file_names.dip_image_path_file],'path_file');
%     end
% else
%     %There is information from a previous run, we know where dip_image is.
%     load([pwd filesep code_names.utilities filesep file_names.dip_image_path_file]);
%     run(path_file);
% end
% clear path_file;
% 
