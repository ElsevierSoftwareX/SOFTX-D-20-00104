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

function TF = calculateTractionForceField(coordinates,faces,stress)

%Keep just the 3 first columns (the 4th one is a label)
faces = faces(:,1:3);

%Reorient the faces so that they all face inwards
faces_reorient = unifyMeshNormals(faces,coordinates,'alignTo','out');

%Compute the normals to the faces at each coordinate
face_normals = STLVertexNormals(faces_reorient, coordinates);
face_normals(isnan(face_normals)) = 0;
n_x = face_normals(:,1);
n_y = face_normals(:,2);
n_z = face_normals(:,3);
clear face_normals;

%Visualization of the normals
% tmp.vertices = coordinates;
% tmp.faces = faces_reorient;
% figure
% patch(tmp,'FaceColor','g','FaceAlpha',0.2), hold on, quiver3(coordinates(:,1),coordinates(:,2),coordinates(:,3),n_x,n_y,n_z), view(3), axis image
% title('Surface normal orientations');

% Compute the number of nodes
n_nodes = 1:size(coordinates,1);

% Compute the traction forces
TF(3*n_nodes-2) = stress(:,1).*n_x + stress(:,4).*n_y + stress(:,5).*n_z;
TF(3*n_nodes-1) = stress(:,4).*n_x + stress(:,2).*n_y + stress(:,6).*n_z;
TF(3*n_nodes) = stress(:,5).*n_x + stress(:,6).*n_y + stress(:,3).*n_z;



% TF(3*numnode-2) = n_x;
% TF(3*numnode-1) = n_y;
% TF(3*numnode) = n_z;