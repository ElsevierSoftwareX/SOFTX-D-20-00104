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

%% Path names

%Name of the path where TFMLAB is stored
path_names.installation_path = pwd;

%Path where elastix is stored
path_names.elastixLibFolder = [path_names.installation_path filesep 'external' filesep 'elastixLib'];

%Or (if your computer is new and IT has blocked elastix)
% path_names.elastixLibFolder = ['C:\Workdir\Programs' filesep 'elastixLib'];

%Path where the templates are stored
path_names.templateShiftCorrectionFilePath = [path_names.installation_path filesep 'cellMat' filesep 'imProcessing' filesep 'shiftCorrection' filesep 'elastix'];
path_names.templateDispCalcFilePath = [path_names.installation_path filesep 'cellMat' filesep 'dispCalc' filesep 'elastix'];

%% Folder names

%Names of the channels
folder_names.cell = 'cell';
folder_names.cellSeg = 'cellSeg';
folder_names.beads = 'beads';
folder_names.fibers = 'fibers';
folder_names.nuclei = 'nuclei';
folder_names.nucleiSeg = 'nucleiSeg';

%Proto-name of the folder in which the raw images will be exported as tiff
folder_names.tiff_raw = 'A_rawImagesTiff';

%Proto-name of the folder in which the filtered images will be stored as
%tiff
folder_names.tiff_filtered = 'B_filteredImages';

%Name of the folder that will contain the shift corrected images
folder_names.shift_corrected = 'C_shiftCorrectedImages';

%Name of the folder where the displacements are stored
folder_names.displacements = 'D_2_displacements';

%Name of the folder where the elastix results are stored
folder_names.elastix_calc = 'D_1_calcElastix';

%Name of the folder where the post-shift corrected images are stored
folder_names.post_shift_corrected = 'D_3_postShiftCorrectedImages';

%Name of the folder where the Paraview files are stored
folder_names.results_for_paraview = 'E_resultsParaview';

%Name of the folder where the tractions are stored
folder_names.tractions = 'F_tractions';

%Name of the subfolders where the traction results are stored
folder_names.abaqus_results_folder = 'results_abaqus';
folder_names.raw_results_folder = 'results_raw';
folder_names.gid_results_folder = 'results_gid';
folder_names.paraview_results_folder = 'results_paraview';
folder_names.intermediate_results_folder = 'results_intermediate';
folder_names.freefem_results_folder = 'results_freefem';
%% File names

%Name of the .mat that is generated to store the filter parameters
file_names.input_param = 'input_param.mat';

%Name of the .mat that is generated to store the filter parameters
file_names.filt_param = 'filt_param.mat';

%Name of the .mat that is generated to store the shift correction
%parameters
file_names.shift_param = 'shift_param.mat';

%Name of the .mat that is generated to store the displacement calculation
%parameters
file_names.disp_param = 'disp_param.mat';

%Name of the .mat that is generated to store the Paraview parameters
file_names.paraview_param = 'paraview_param.mat';

%Name of the .mat that is generated to store the traction calculation
%parameters
file_names.trac_param = 'trac_param.mat';

%Name of the paraview template
file_names.paraview_template = 'template.pvsm';

%Name of the file that contains the dip_image path
file_names.dip_image_path_file = 'dip_image_path.mat';

%% Shift correction parameters

% Parameters for 'elastix'
shiftParam.elastix.multiscale.num = 3;
shiftParam.elastix.multiscale.pyramidSchedule = [4 4 4;  2 2 2;  1 1 1]; % For each scale; [x,y,z]
shiftParam.elastix.metric.type = 'mutualInfo'; % 'normXCorr' or 'mutualInfo' or 'meanSquares'
shiftParam.elastix.metric.mutualInfo.numHistBins = [16 32 32]; % For each scale. Typically 16 or 32 are good values. Only valid if (shiftParam.pre.elastix.metric.type='mutualInfo')
shiftParam.elastix.optim.method = 'adapStochGradDesc'; % 'qnLBFGS' or 'adapStochGradDesc'
shiftParam.elastix.optim.asgd.speed = false; % Only valid if (shiftParam.pre.elastix.optim.method='adapStochGradDesc')
shiftParam.elastix.optim.iterNum = [20 40 80];% For each scale.
shiftParam.elastix.optim.evalSampleNum = [5000 5000 5000]; % For each scale. Empty [] --> full 

%% FFD parameters
dispParam.grid.minSpacing.xy = 5; % units in um (positive integer number )- average distance between two beads approx
dispParam.grid.minSpacing.z = 5; % units in um (only needed for 3D Analysis -- positive integer number)
dispParam.mask.beads = false; % true or false, leave as false
dispParam.mask.cell = true; % true or false - always leave true, ignores the beads inside the mask of the cell
dispParam.strategy = 'direct'; % 'direct' or 'sequential' or 'seqBackProp' or 'bwdIniTf'
dispParam.artifactCorr = false; % true, false

% Advanced parameters:
dispParam.refTimePoint = []; % which timepoint will be used as relaxed state (if empty [], the last timepoint will be used) --> only valid if dispParam.strategy ='direct' and a file for the relaxed state has not been provided
dispParam.seq.tStep =1; % Only valid if (dispParam.strategy = 'seqBackProp' or 'sequential') you can define timestep for comparing sequential data and choosing reference state
dispParam.multiscale.num = 3;
dispParam.multiscale.pyramidSchedule = [4 4 4;  2 2 2;  1 1 1]; % For each scale; [x,y,z]
dispParam.grid.spacingSchedule = [4 4 4; 2 2 2; 1 1 1]; % For each scale; [x,y,z]
dispParam.grid.passiveEdgeGridPointNum = [0 0 0] ; % For each scale <-- An integer value in the range [0,4]. If not sure, use 0. A value of 4 will avoid all deformations at the edge of the image. Make sure that 2*passiveEdgeGridPointNum < ControlPointGridSize in each dimension
dispParam.metric.type = 'normXCorr'; % 'normXCorr' or 'mutualInfo' or 'meanSquares'
dispParam.metric.mutualInfo.numHistBins = [16 32 32]; % For each scale. Typically 16 or 32 are good values. Only valid if (shiftParam.pre.elastix.metric.type='mutualInfo')
dispParam.optim.method = 'adapStochGradDesc'; % 'qnLBFGS' or 'adapStochGradDesc'
dispParam.optim.asgd.speed = false; % Only valid if (dispParam.optim.method='adapStochGradDesc')
dispParam.optim.iterNum = [500 500 500];% For each scale.
dispParam.optim.evalSampleNum = [5000 5000 5000]; % For each scale. Empty [] --> full 

%% Code names
code_names.guis = 'guis';
code_names.images = 'images';
code_names.external = 'external';
code_names.utilities = 'utilities';
code_names.microscope2Tif = 'microscope2Tif';
code_names.cellMat = 'cellMat';
code_names.results2Paraview = 'results2Paraview';
code_names.tractionRecovery = 'tractionRecovery';


