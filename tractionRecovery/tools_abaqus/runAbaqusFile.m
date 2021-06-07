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
function runAbaqusFile(write_options,mech_props,mode)

%Save the current path
code_path = pwd;

%First go to where the file is
cd(write_options.write_path);

%Second, run Abaqus file
tmp1 = ['"' write_options.abq1_text '" amd64'];
tmp2 = ['"' write_options.abq2_text '" intel64'];
switch mode
    case 'stiffness_matrix'
        switch mech_props.ecm_behavior
            case 'Non linear elastic'
                command_line = [tmp1 ' & ' tmp2 ' & abaqus j=',write_options.write_filename,' user=wrt_rf.for'];
            case 'Linear elastic'
                command_line = ['abaqus j=' write_options.write_filename];
        end
    case 'results'
        %command_line = ['abq1 & abq2 & abaqus j=',write_options.write_filename,'_check user=wrt_us.for'];
        command_line = [tmp1 ' & ' tmp2 ' & abaqus j=',write_options.write_filename,'_check user=wrt_us.for'];
end
clear tmp1 tmp2;
status = system(command_line);
pause(4);

% Wait for file to be created.
maxSecondsToWait = 36000; % Wait 10 hours at most.
secondsWaitedSoFar  = 0;
while secondsWaitedSoFar < maxSecondsToWait    
    switch mode
        case 'stiffness_matrix'
            condition = ~exist([write_options.write_filename '.lck'], 'file');
        case 'results'
            condition = ~exist([write_options.write_filename,'_check.lck'], 'file');
    end        
    if condition
        break;
    end       
    pause(1); % Wait 1 second.
    secondsWaitedSoFar = secondsWaitedSoFar + 1;     
end

%Go back to the code path
cd(code_path);