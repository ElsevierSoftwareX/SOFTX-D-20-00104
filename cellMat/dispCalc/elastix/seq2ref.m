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

function seq2ref(timePoints,seqTstep,calcFolder,elastixLibFolder,tformFolderBaseName,tformRunTemplate,matFileBaseName,resolution,outImFormat)
% This function converts all the displacements that have been calculated
% relatively (between consecutive time points) to "absolute" displacements
% that has the relaxed state as the reference


tpNum = length(timePoints);



%% Build the required folders

% Main folder for concatenations
concFolder = [calcFolder filesep 'tformConcatenation'];
if ~exist(concFolder,'dir')
    mkdir(concFolder)
end
% Folder for the images
imFolder = [concFolder filesep 'images'];
if ~exist(imFolder,'dir')
    mkdir(imFolder)
end
% Folder for each timepoint
for tt=1:tpNum
    if ~exist([concFolder filesep 't' num2str(timePoints(tt))],'dir')
        mkdir([concFolder filesep 't' num2str(timePoints(tt))])
    end
end



%% Process the timelapse


movingImFileName = 'movingIm';
imFormat = 'tiff'; % default format


for tt=1:tpNum
    
    tp = timePoints(tt);
    disp(' ... ... ...')
    disp(['Backward Prop. Transform Concatenation : TIME POINT ' num2str(tp) ' ... ...'])
    
    currentFolder = [concFolder filesep 't' num2str(timePoints(tt))];
    
    % Generate the concatenated transforms
    for jj=tt:seqTstep:tpNum
        sourceFile = [tformFolderBaseName '_t' num2str(timePoints(jj)) filesep 'TransformParameters.0.txt'];
        destinationFile = [currentFolder filesep 'conc' num2str(timePoints(jj)) 'TransformParameters.0.txt'];
        if jj==tt
            copyfile(sourceFile,destinationFile);
        else
            strIn{1} = 'NoInitialTransform'; strOut{1} = [currentFolder filesep 'conc' num2str(timePoints(jj-seqTstep)) 'TransformParameters.0.txt'];
            replaceMultiString(strIn,strOut,sourceFile,destinationFile)
            clear strIn strOut
        end
    end
    clear sourceFile destinationFile
    
    % Copy/Modify the running (command line) file
    currentRunFile = [currentFolder filesep 'elastixConcTransformRun.bat'];
    strIn{1} = 'elastixLibPath'; strOut{1} = elastixLibFolder;
    strIn{2} = 'imReadPath'; strOut{2} = imFolder; 
    strIn{3} = 'tformReadPath'; strOut{3} = currentFolder; 
    strIn{4} = 'savePath'; strOut{4} = currentFolder;  
    strIn{5} = 'movingImName'; strOut{5} = [movingImFileName '.' imFormat];
    strIn{6} = 'tformFileName'; strOut{6} = ['conc' num2str(timePoints(max(tt:seqTstep:tpNum))) 'TransformParameters.0.txt'];
    replaceMultiString(strIn,strOut,tformRunTemplate,currentRunFile);
    clear strIn strOut  
    
    % Write the image to a file
    load([matFileBaseName '_t' num2str(tp) '.mat'],'stressedIm');
    if exist([imFolder filesep movingImFileName '.' imFormat],'file')
        delete([imFolder filesep movingImFileName '.' imFormat])
    end
    writeIm2File(stressedIm,imFolder,movingImFileName,imFormat);
    clear stressedIm
    
    % Calculate the displacements
    [dispField,registeredIm] = ffdImReg(currentRunFile,currentFolder,outImFormat,resolution,['conc' num2str(timePoints(max(tt:seqTstep:tpNum))) 'TransformParameters.0.txt']);
    % Store the results
    save([matFileBaseName '_t' num2str(tp) '.mat'],'-append','dispField','registeredIm')
    
end


% Delete some of the folders
try
    rmdir(imFolder,'s')
catch
    disp('Impossible to delete the temporal folder containing moving images during the concatenation of the transforms')
end
