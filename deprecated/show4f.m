function show(elements,coordinates,uf,filename,bf,tf,ssf,srf)
%SHOW  Plots three-dimensional solution
%    SHOW(ELEMENTS,DIRICHLET,NEUMANN,COORDINATES,U,LAMBDA,MU) plots the
%    strained mesh and visualizes the stresses in grey tones.
%
%
%   See also FEM_LAME3D.

%    J. Alberty, C. Carstensen, S. A. Funken, and R. Klose  07-03-00
%    File <show.m> in $(HOME)/acfk/fem_lame3d/angle/

% E=zeros(4*size(elements,1),3);
%  C=zeros(size(elements,1),1);
%  for j=1:size(elements,1)
%    PhiGrad=[1,1,1,1;coordinates(elements(j,:),:)']\[zeros(1,3);eye(3)]; 
%    U_Grad = u([1;1;1]*3*elements(j,:)-[2;1;0]*[1,1,1,1])*PhiGrad;
%    SIGMA = lambda * trace(U_Grad)*eye(3)+mu*(U_Grad + U_Grad');
%    C(j) = sqrt(sum(eig(SIGMA).^2));   
%  end
%  Area = zeros(size(elements,1),1);
%  AreaOmega = zeros(max(max(elements)),1);
%  AvC = zeros(max(max(elements)),1);
%  for j=1:size(elements,1)
%    Area = det([1,1,1,1;coordinates(elements(j,:),:)'])/6;
%    AreaOmega(elements(j,:))=AreaOmega(elements(j,:))+Area;
%    AvC(elements(j,:)) = AvC(elements(j,:))+Area*[1;1;1;1]*C(j);
%  end
%  AvC = AvC./AreaOmega;
%  E=[dirichlet;neumann];
%   factor=100.0;
%   colormap(1-gray);
%   trisurf(E,factor*u(1:3:size(u,1))+coordinates(:,1), ...
%         factor*u(2:3:size(u,1))+coordinates(:,2), ...
%         factor*u(3:3:size(u,1))+coordinates(:,3),AvC,'facecolor','interp');
%   view(-50,35)
%   axis equal; axis([-10 70 -100 20 0 40])
%   xlabel('x'); ylabel('y'); zlabel('z');
  



% factor=0;
% %colormap(1-gray)
% trisurf(elements,factor*u(1:3:size(u,1))+coordinates(:,1), ...
%     factor*u(2:3:size(u,1))+coordinates(:,2),factor*u(3:3:size(u,1))+coordinates(:,3), ...
%     u(1:2:size(u,1)), 'facecolor','interp');
% view(0,90)
% colorbar('vert')






vc = 1:size(coordinates,1);

fid = fopen(strcat(filename,'.flavia.res'),'wt');
fprintf(fid,'Gid Post Results File 1.0');
fprintf(fid,'\n','');
fprintf(fid,'RESULT "Displacements Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[vc', uf(1:3:size(uf,1)), uf(2:3:size(uf,1)), uf(3:3:size(uf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');


fprintf(fid,'RESULT "Volume Forces Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[vc', bf(1:3:size(bf,1)), bf(2:3:size(bf,1)), bf(3:3:size(bf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');


fprintf(fid,'RESULT "Traction Forces Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[vc', tf(1:3:size(tf,1)), tf(2:3:size(tf,1)), tf(3:3:size(tf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

fprintf(fid,'RESULT "Principal Stress Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[vc', ssf(:,1), ssf(:,2), ssf(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');


fprintf(fid,'RESULT "Principal Strain Forward" "DD" 1 Vector OnNodes');
fprintf(fid,'\n','');
fprintf(fid,'VALUES');
fprintf(fid,'\n','');
rd=[vc', srf(:,1), srf(:,2), srf(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fprintf(fid,'END VALUES');
fprintf(fid,'\n','');

fclose(fid);




