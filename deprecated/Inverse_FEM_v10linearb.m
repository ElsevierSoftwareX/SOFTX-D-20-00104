clear 
close all
clc




for i=7:9

filename = strcat('FI8_FB_t',num2str(i),'_void')    
fln2 = strcat('F8_t',num2str(i),'.mat')


load(fln2);

%%%%%%%%%%%%% Carga de datos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('F8_morphoDispMetrics.mat');

mask_all = cellMask.mask3D;





Lx = param.res.xy;
Ly = param.res.xy;
Lz = param.res.z;       


for k=1:size(mask_all,3)
    for j=1:round(Lz/Lx,0)
        mask_all2(:,:,round(Lz/Lx,0)*(k-1)+j) = mask_all(:,:,k);
    end
end

%%%%%%%%%%%%%%% Paso a malla de tetraedros iso2mesh %%%%%%%%%%%%%%%%%%%%%%%


%addpath('C:\Program Files\MATLAB\R2017a\toolbox\iso2mesh');
addpath('C:\Users\r0726821\KUL_jsanz\iso2mesh');

[no2,el2,regions,holes]=v2s(mask_all2,0.5,3);

%addpath('C:\Users\Erreina\Downloads\smoothpatch_version1b');
addpath('C:\Users\r0726821\KUL_jsanz\smoothpatch_version1b');
AA.vertices=no2;
AA.faces=el2(:,1:3);
AA2=smoothpatch(AA,1,5);





IID = find(mask_all2==1);

[ID JD KD] = ind2sub(size(mask_all2),IID);

%[node,el1,fc1]=surf2mesh(AA2.vertices,AA2.faces,[0 -1 -1],[size(mask_all2,1) size(mask_all2,2)+1 size(mask_all2,3)+1],1,100,[ID(1) JD(1) KD(1)],[],1);
%[node,el1,fc1]=surf2mesh(AA2.vertices,AA2.faces,[17 -1 -20],[size(mask_all2,1)+37 size(mask_all2,2) size(mask_all2,3)+20],1,100,[ID(1) JD(1) KD(1)],[],1);
[node,el1,fc1]=surf2mesh(AA2.vertices,AA2.faces,[-1 -1 -1],[size(mask_all2,1)+1 size(mask_all2,2)+1 size(mask_all2,3)+1],1,200,[ID(1) JD(1) KD(1)],[],1);




elem0=el1((el1(:,5)==0),:);
elem=el1((el1(:,5)==1),:);

ind_e0 = find(el1(:,5)==0);
ind_e1 = find(el1(:,5)==1);

node = [Lx*node(:,1) Ly*node(:,2) Lz/round(Lz/Lx,0)*node(:,3)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Fichero para GID

fid = fopen(strcat(filename,'.flavia.msh'),'wt');
fprintf(fid,'MESH  "Sprout" DIMENSION 3 ELEMTYPE Tetrahedra Nnode 4');
fprintf(fid,'\n','');
fprintf(fid,'COORDINATES');
fprintf(fid,'\n','');


rd=[(1:size(node,1))', node];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');

fprintf(fid,'END COORDINATES');
fprintf(fid,'\n','');
fprintf(fid,'ELEMENTS');
fprintf(fid,'\n','');



rr=[(size(elem0,1)+1:size(elem0,1)+size(elem,1))', elem(:,1), elem(:,2), elem(:,3), elem(:,4)];
fprintf(fid,'%7i %7i %7i %7i %7i\n',rr');

	
fprintf(fid,'END ELEMENTS');
fprintf(fid,'\n','');

fprintf(fid,'MESH  "Gel" DIMENSION 3 ELEMTYPE Tetrahedra Nnode 4');
fprintf(fid,'\n','');
fprintf(fid,'COORDINATES');
fprintf(fid,'\n','');
fprintf(fid,'END COORDINATES');
fprintf(fid,'\n','');
fprintf(fid,'ELEMENTS');
fprintf(fid,'\n','');

rr=[(1:size(elem0,1))', elem0(:,1), elem0(:,2), elem0(:,3), elem0(:,4)];
fprintf(fid,'%7i %7i %7i %7i %7i\n',rr');


fprintf(fid,'END ELEMENTS');
fprintf(fid,'\n','');

fclose(fid);


%%%%%%%%%%%%%%%%% Campo de desplazamientos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Nx Ny Nz] = size(dispField.Y);

u_read = reshape(dispField.Y,Nx*Ny*Nz,1);
v_read = reshape(dispField.X,Nx*Ny*Nz,1);
w_read = reshape(dispField.Z,Nx*Ny*Nz,1);


[ix,iy,iz] = ind2sub([Nx Ny Nz],1:(Nx)*(Ny)*(Nz));

x_read(sub2ind([Nx Ny Nz],ix,iy,iz)) = ix*Lx;
y_read(sub2ind([Nx Ny Nz],ix,iy,iz)) = iy*Ly;
z_read(sub2ind([Nx Ny Nz],ix,iy,iz)) = iz*Lz;

% x_read = x_read+(min(node(:,1))-min(x_read))+Lx/2;
% y_read = y_read+(min(node(:,2))-min(y_read))+Ly/2;
% z_read = z_read+(min(node(:,3))-min(z_read))+Ly/2;

FUX = scatteredInterpolant(x_read',y_read',z_read',double(u_read));
UX = FUX(node(:,1),node(:,2),node(:,3));
clear FUX

FUY = scatteredInterpolant(x_read',y_read',z_read',double(v_read));
UY = FUY(node(:,1),node(:,2),node(:,3));
clear FUY

FUZ = scatteredInterpolant(x_read',y_read',z_read',double(w_read));
UZ = FUZ(node(:,1),node(:,2),node(:,3));
clear FUZ

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coordinates = [node(:,1),node(:,2),node(:,3)];

u(1:3:3*size(coordinates,1)) = UX;
u(2:3:3*size(coordinates,1)) = UY;
u(3:3:3*size(coordinates,1)) = UZ;

elements = [elem0;elem];




%%%%%%%%%%%%%%%%%%%%%%% Output Abaqus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(strcat(filename,'.inp'),'wt');
fprintf(fid,'*HEADING');
fprintf(fid,'\n','');
fprintf(fid,'*NODE,NSET=NTODOS');
fprintf(fid,'\n','');

nnode = 1:size(coordinates,1);

rd=[nnode', coordinates(:,1), coordinates(:,2), coordinates(:,3)];
fprintf(fid,'%7i, %12.7f, %12.7f, %12.7f\n',rd');


fprintf(fid,'*ELEMENT,TYPE=C3D4,ELSET=GEL');
fprintf(fid,'\n','');

rr=[(1:size(elem0,1))', elem0(:,1), elem0(:,2), elem0(:,3), elem0(:,4)];
fprintf(fid,'%7i, %7i, %7i, %7i, %7i\n',rr');

fprintf(fid,'*ELEMENT,TYPE=C3D4,ELSET=SPROUT');
fprintf(fid,'\n','');

clear rr

rr=[(size(elem0,1)+1:size(elem0,1)+size(elem,1))', elem(:,1), elem(:,2), elem(:,3), elem(:,4)];
fprintf(fid,'%7i, %7i, %7i, %7i, %7i\n',rr');


fprintf(fid,'*SOLID SECTION, ELSET=SPROUT, MATERIAL=MATPROPS2');
fprintf(fid,'\n','');
fprintf(fid,'*SOLID SECTION, ELSET=GEL, MATERIAL=MATPROPS');
fprintf(fid,'\n','');
fprintf(fid,'*MATERIAL,NAME=MATPROPS2'); 
fprintf(fid,'\n','');
fprintf(fid,'*ELASTIC,TYPE=ISOTROPIC');
fprintf(fid,'\n','');
fprintf(fid,'100.0e-6,  4.50000E-01');
fprintf(fid,'\n','');
fprintf(fid,'*MATERIAL,NAME=MATPROPS'); 
fprintf(fid,'\n','');
fprintf(fid,'*ELASTIC,TYPE=ISOTROPIC');
fprintf(fid,'\n','');
fprintf(fid,'100.0,  4.50000E-01');
fprintf(fid,'\n','');
fprintf(fid,'*STEP'); 
fprintf(fid,'\n','');
fprintf(fid,'*MATRIX GENERATE, STIFFNESS');
fprintf(fid,'\n','');
fprintf(fid,'*MATRIX OUTPUT, STIFFNESS, FORMAT=COORDINATE');
fprintf(fid,'\n','');
fprintf(fid,'*END STEP');
fprintf(fid,'\n','');


fclose(fid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


command2 = strcat('abaqus j=',filename);
status=dos(command2)


 pause(10);

% Wait for file to be created.
maxSecondsToWait = 36000; % Wait 10 hours at most.
secondsWaitedSoFar  = 0;
while secondsWaitedSoFar < maxSecondsToWait 
  if ~exist(strcat(filename,'.lck'), 'file')
    break;
  end
  pause(1); % Wait 1 second.
  secondsWaitedSoFar = secondsWaitedSoFar + 1;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asparse=load(strcat(filename,'_STIF1.mtx'));
A = spconvert(asparse);



% RESOLVER

t0 = A*u';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%% Problema inverso %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ind_fc = find(fc1(:,4)==0);
% 
% set1 = unique(fc1(ind_fc,1:3));
% 
% neig_1 = vecino(set1,elem0(:,1:4));
% neig_2 = vecino(neig_1,elem0(:,1:4));
% neig_3 = vecino(neig_2,elem0(:,1:4));

% ind = unique([unique(elem(:,1:4));unique(fc1(:,1:3));neig_2]);

ind = unique([unique(elem(:,1:4));unique(fc1(:,1:3))]);


ind_a = [3*ind-2;3*ind-1;3*ind];
ind_b = setdiff(1:length(u),ind_a);

ua_s = u(ind_a)';
ub_s = u(ind_b)';

ia = length(ind_a);
ib = length(ind_b);

% x = gmres(@(x)evaluador3(x,A,ind_a,ind_b),[ua_s; ub_s; zeros(length(ind_b),1)],[],1e-6,3000);

%%%%%%%%%%%%% funciona y converge - tarda 5 horas y media 
x = gmres(@(x)evaluador3(x,A,ind_a,ind_b),[ua_s; ub_s; zeros(length(ind_b),1)],[],1e-3,30000);

% x = pcg(@(x)evaluador3(x,A,ind_a,ind_b),[ua_s; ub_s; zeros(length(ind_b),1)],1e-6,30000);

% AA = [speye(ia) speye(ia,ib)-speye(ia,ib) A(ind_a,ind_b); speye(ib,ia)-speye(ib,ia) speye(ib) A(ind_b,ind_b); A(ind_b,ind_a) A(ind_b,ind_b) speye(ib,ib)-speye(ib,ib)];
% % [L,U,P] = lu(AA);
% [L,U] = ilu(AA,struct('type','ilutp','droptol',1e-3));
% x = gmres(@(x)evaluador3(x,A,ind_a,ind_b),[ua_s; ub_s; zeros(length(ind_b),1)],[],1e-3,3000,L,U);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/UMFPACK/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/CHOLMOD/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/AMD/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/COLAMD/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/CCOLAMD/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/CAMD/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/ssget
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/CXSparse/MATLAB/Demo
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/CXSparse/MATLAB/CSparse
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/LDL/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/BTF/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/KLU/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/SPQR/MATLAB
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/RBio/RBio
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/Factorize
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/MESHND
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/LINFACTOR
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/find_components
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/GEE
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/shellgui
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/waitmex
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/spqr_rank
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/spqr_rank/SJget
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/UFcollection
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/SSMULT
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/dimacs10
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/spok
% addpath C:\Users\r0726821\KUL_jsanz\Sprout\SuiteSparse-5.2.0\SuiteSparse/MATLAB_Tools/sparseinv
% 
% bb = [ua_s; ub_s; zeros(length(ind_b),1)];
% AA = [speye(ia) speye(ia,ib)-speye(ia,ib) A(ind_a,ind_b); speye(ib,ia)-speye(ib,ia) speye(ib) A(ind_b,ind_b); A(ind_b,ind_a) A(ind_b,ind_b) speye(ib,ib)-speye(ib,ib)];
% x = spqr_solve(AA,bb);






% bb = [ua_s; ub_s; zeros(length(ind_b),1)];
% AA = [speye(ia) speye(ia,ib)-speye(ia,ib) A(ind_a,ind_b); speye(ib,ia)-speye(ib,ia) speye(ib) A(ind_b,ind_b); A(ind_b,ind_a) A(ind_b,ind_b) speye(ib,ib)-speye(ib,ib)];
% x = AA\bb;


uu(ind_a) = x(1:length(ind_a));
uu(ind_b) = x(length(ind_a)+1:length(ind_a)+length(ind_b));

tt = A*uu';


%%%%%%%%%%%%%%%%%%%%%%%% Control Abaqus inverse %%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(strcat(filename,'_check.inp'),'wt');
fprintf(fid,'*HEADING');
fprintf(fid,'\n','');
fprintf(fid,'*NODE,NSET=NTODOS');
fprintf(fid,'\n','');

nnode = 1:size(coordinates,1);

rd=[nnode', coordinates(:,1), coordinates(:,2), coordinates(:,3)];
fprintf(fid,'%7i, %12.7f, %12.7f, %12.7f\n',rd');


fprintf(fid,'*ELEMENT,TYPE=C3D4,ELSET=GEL');
fprintf(fid,'\n','');

rr=[(1:size(elem0,1))', elem0(:,1), elem0(:,2), elem0(:,3), elem0(:,4)];
fprintf(fid,'%7i, %7i, %7i, %7i, %7i\n',rr');

fprintf(fid,'*ELEMENT,TYPE=C3D4,ELSET=SPROUT');
fprintf(fid,'\n','');

clear rr

rr=[(size(elem0,1)+1:size(elem0,1)+size(elem,1))', elem(:,1), elem(:,2), elem(:,3), elem(:,4)];
fprintf(fid,'%7i, %7i, %7i, %7i, %7i\n',rr');


fprintf(fid,'*SOLID SECTION, ELSET=SPROUT, MATERIAL=MATPROPS2');
fprintf(fid,'\n','');
fprintf(fid,'*SOLID SECTION, ELSET=GEL, MATERIAL=MATPROPS');
fprintf(fid,'\n','');
fprintf(fid,'*MATERIAL,NAME=MATPROPS2'); 
fprintf(fid,'\n','');
fprintf(fid,'*ELASTIC,TYPE=ISOTROPIC');
fprintf(fid,'\n','');
fprintf(fid,'100.0e-6,  4.50000E-01');
fprintf(fid,'\n','');
fprintf(fid,'*MATERIAL,NAME=MATPROPS'); 
fprintf(fid,'\n','');
fprintf(fid,'*ELASTIC,TYPE=ISOTROPIC');
fprintf(fid,'\n','');
fprintf(fid,'100.0,  4.50000E-01');
fprintf(fid,'\n','');
fprintf(fid,'*STEP,INC=1000');
fprintf(fid,'\n','');
fprintf(fid,'*STATIC,DIRECT');
fprintf(fid,'\n','');
fprintf(fid,'1.0,1.0');
fprintf(fid,'\n','');
fprintf(fid,'*BOUNDARY,TYPE=DISPLACEMENT');
fprintf(fid,'\n','');

clear rr

indi = nnode;

abb = 1:size(coordinates,1);

rr=[abb',ones(length(indi),2),u(3*indi-2)'; abb',2*ones(length(indi),2),u(3*indi-1)'; abb', 3*ones(length(indi),2), u(3*indi)'];
fprintf(fid,'%7i, %7i, %7i, %12.7f\n',rr');
fprintf(fid,'*NODE FILE');
fprintf(fid,'\n','');
fprintf(fid,'U');
fprintf(fid,'\n','');
fprintf(fid,'*EL FILE,ELSET=GEL,POSITION=AVERAGED AT NODES');
fprintf(fid,'\n','');
fprintf(fid,'S,SP,EP');
fprintf(fid,'\n','');
fprintf(fid,'*END STEP');
fprintf(fid,'\n','');


fprintf(fid,'*STEP,INC=1000');
fprintf(fid,'\n','');
fprintf(fid,'*STATIC,DIRECT');
fprintf(fid,'\n','');
fprintf(fid,'1.0,1.0');
fprintf(fid,'\n','');
fprintf(fid,'*BOUNDARY,TYPE=DISPLACEMENT');
fprintf(fid,'\n','');

clear rr

indi = nnode;

abb = 1:size(coordinates,1);

rr=[abb',ones(length(indi),2),uu(3*indi-2)'; abb',2*ones(length(indi),2),uu(3*indi-1)'; abb', 3*ones(length(indi),2), uu(3*indi)'];
fprintf(fid,'%7i, %7i, %7i, %12.7f\n',rr');
fprintf(fid,'*NODE FILE');
fprintf(fid,'\n','');
fprintf(fid,'U');
fprintf(fid,'\n','');
fprintf(fid,'*EL FILE,ELSET=GEL,POSITION=AVERAGED AT NODES');
fprintf(fid,'\n','');
fprintf(fid,'S,SP,EP');
fprintf(fid,'\n','');
fprintf(fid,'*END STEP');
fprintf(fid,'\n','');


fclose(fid)


command2 = strcat('abaqus j=',filename,'_check user=wrt_us.for');
status=dos(command2)


 pause(10);

% Wait for file to be created.
maxSecondsToWait = 36000; % Wait 10 hours at most.
secondsWaitedSoFar  = 0;
while secondsWaitedSoFar < maxSecondsToWait 
  if ~exist(strcat(filename,'_check.lck'), 'file')
    break;
  end
  pause(1); % Wait 1 second.
  secondsWaitedSoFar = secondsWaitedSoFar + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%% Visualizacion GID %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BDF_i = volumeforces(coordinates,elements,tt);

BDF_f = volumeforces(coordinates,elements,t0);


stress_f(1:size(node,1),1:6) = 0;
stress_i(1:size(node,1),1:6) = 0;
stressp_f(1:size(node,1),1:3) = 0;
stressp_i(1:size(node,1),1:3) = 0;
strainp_f(1:size(node,1),1:3) = 0;
strainp_i(1:size(node,1),1:3) = 0;


load('stresses.dat');
load('stressesp.dat');
load('strainsp.dat');

stress_f(unique(elem0(:,1:4)),:) = stresses(1:length(unique(elem0(:,1:4))),2:end);
stress_i(unique(elem0(:,1:4)),:) = stresses(length(unique(elem0(:,1:4)))+1:end,2:end);

stressp_f(unique(elem0(:,1:4)),:) = stressesp(1:length(unique(elem0(:,1:4))),2:end);
stressp_i(unique(elem0(:,1:4)),:) = stressesp(length(unique(elem0(:,1:4)))+1:end,2:end);

strainp_f(unique(elem0(:,1:4)),:) = strainsp(1:length(unique(elem0(:,1:4))),2:end);
strainp_i(unique(elem0(:,1:4)),:) = strainsp(length(unique(elem0(:,1:4)))+1:end,2:end);


TF_i = tractionforces2(coordinates,elements,fc1,stress_i,elem0);
TF_f = tractionforces2(coordinates,elements,fc1,stress_f,elem0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vc = 1:size(coordinates,1);

uf = u';
fid = fopen(strcat('uf',fln2),'wt');
rd=[vc', uf(1:3:size(uf,1)), uf(2:3:size(uf,1)), uf(3:3:size(uf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)

ui = uu';
fid = fopen(strcat('ui',fln2),'wt');
rd=[vc', ui(1:3:size(uf,1)), ui(2:3:size(uf,1)), ui(3:3:size(uf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)

bf =BDF_f';
fid = fopen(strcat('bf',fln2),'wt');
rd=[vc', bf(1:3:size(bf,1)), bf(2:3:size(bf,1)), bf(3:3:size(bf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)

bi =BDF_i';
fid = fopen(strcat('bi',fln2),'wt');
rd=[vc', bi(1:3:size(bf,1)), bi(2:3:size(bf,1)), bi(3:3:size(bf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)

tf =TF_f';
fid = fopen(strcat('tf',fln2),'wt');
rd=[vc', tf(1:3:size(tf,1)), tf(2:3:size(tf,1)), tf(3:3:size(tf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)

ti =TF_i';
fid = fopen(strcat('ti',fln2),'wt');
rd=[vc', ti(1:3:size(tf,1)), ti(2:3:size(tf,1)), ti(3:3:size(tf,1))];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)

ssf =stressp_f;
fid = fopen(strcat('ssf',fln2),'wt');
rd=[vc', ssf(:,1), ssf(:,2), ssf(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)

ssi =stressp_i;
fid = fopen(strcat('ssi',fln2),'wt');
rd=[vc', ssi(:,1), ssi(:,2), ssi(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)

srf =strainp_f;
fid = fopen(strcat('srf',fln2),'wt');
rd=[vc', srf(:,1), srf(:,2), srf(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)

sri =strainp_i;
fid = fopen(strcat('sri',fln2),'wt');
rd=[vc', sri(:,1), sri(:,2), sri(:,3)];
fprintf(fid,'%7i %12.7f %12.7f %12.7f\n',rd');
fclose(fid)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

show4(elements,coordinates,u',uu',filename,BDF_f',BDF_i',TF_f',TF_i',stressp_f,strainp_f,stressp_i,strainp_i); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%   Indicadores  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Desplazamientos
% 
% ud = u;
% td = t0;
% 
% ub = load('u_base.mat');
% tb = load('t_base.mat');
% 
% k=0;
% epsi = 1e-9;
% for i=3:3:length(u)
%     k=k+1;
%     ub_m(k) = sqrt(ub.u(i-2)^2+ub.u(i-1)^2+ub.u(i)^2);
%     tb_m(k) = sqrt(tb.t0(i-2)^2+tb.t0(i-1)^2+tb.t0(i)^2);
%     
%     ud_m(k) = sqrt(ud(i-2)^2+ud(i-1)^2+ud(i)^2);
%     td_m(k) = sqrt(td(i-2)^2+td(i-1)^2+td(i)^2);
%     
%     ui_m(k) = sqrt(uu(i-2)^2+uu(i-1)^2+uu(i)^2);
%     ti_m(k) = sqrt(tt(i-2)^2+tt(i-1)^2+tt(i)^2);
% end
% 
% 
% err_ti = 0;
% err_td = 0;
% k = 0;
% for i=1:length(ub_m)
%     if abs(tb_m(i))>epsi
%         err_ti = err_ti + abs((ti_m(i)-tb_m(i))/tb_m(i));
%         err_td = err_td + abs((td_m(i)-tb_m(i))/tb_m(i));
%         k = k+1;
%     end
% end
% 
% err_ti = err_ti*100/k
% err_td = err_td*100/k
% 
% err_ui = 0;
% err_ud = 0;
% k = 0;
% for i=1:length(ub_m)
%     if abs(ub_m(i))>epsi
%         err_ui = err_ui + abs((ui_m(i)-ub_m(i))/ub_m(i));
%         err_ud = err_ud + abs((ud_m(i)-ub_m(i))/ub_m(i));
%         k = k+1;
%     end
% end
% 
% err_ui = err_ui*100/k
% err_ud = err_ud*100/k

end

