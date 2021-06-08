clear 
close all
clc


fln2 = 'int_05';


filename=strcat('FINL_',fln2);
filename2='Prueba_Single_Cell_void';  

%%%%%%%%%%%%% Carga de datos %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(strcat(filename2,'_node.mat'));
load(strcat(filename2,'_el1.mat'));
load(strcat(filename2,'_fc1.mat'));

elem0=el1((el1(:,5)==1),:);
elem=el1((el1(:,5)==0),:);

ind_e0 = find(el1(:,5)==1);
ind_e1 = find(el1(:,5)==0);

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

load('disp_corr.mat');

UX = disp_corr(:,1);
UY = disp_corr(:,2);
UZ = disp_corr(:,3);

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
fprintf(fid,'*Hyperelastic, neo hooke');
fprintf(fid,'\n','');
fprintf(fid,'16.72241,6.0e-04');
fprintf(fid,'\n','');
fprintf(fid,'*STEP,NLGEOM,INC=1000');
fprintf(fid,'\n','');
fprintf(fid,'*STATIC');
fprintf(fid,'\n','');
fprintf(fid,'1.0,1.0,0.025,1.0');
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
fprintf(fid,'RF');
fprintf(fid,'\n','');
fprintf(fid,'*END STEP');
fprintf(fid,'\n','');
fprintf(fid,'*STEP,NLGEOM');
fprintf(fid,'\n','');
fprintf(fid,'*MATRIX GENERATE, STIFFNESS, ELEMENT BY ELEMENT');
fprintf(fid,'\n','');
fprintf(fid,'*MATRIX OUTPUT, STIFFNESS, FORMAT=COORDINATE');
fprintf(fid,'\n','');
fprintf(fid,'*END STEP');


fclose(fid)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


command2 = strcat('abq1 & abq2 & abaqus j=',filename,' user=wrt_rf.for');
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

bsparse=load(strcat(filename,'_STIF2.mtx'));

ind_i = bsparse(:,1)*12 + bsparse(:,2);
ind_j = bsparse(:,3);


BSF = spconvert([ind_i ind_j bsparse(:,4)]);


KS1 = (reshape(BSF',12*12,size(BSF,1)/12))';

clear Ig
clear Jg 

Ig = 3*elements(:,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
for i=1:11
    Ig = [Ig,3*elements(:,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0]];
end

for i=1:4
    if i==1
        Jg = 3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[2,2,2,2,2,2,2,2,2,2,2,2];
    else
        Jg = [Jg,3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[2,2,2,2,2,2,2,2,2,2,2,2]];
    end
    Jg = [Jg,3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[1,1,1,1,1,1,1,1,1,1,1,1]];
    Jg = [Jg,3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[0,0,0,0,0,0,0,0,0,0,0,0]];
end


AA = sparse ( Ig , Jg , KS1, 3*size(coordinates,1) , 3*size(coordinates,1) ) ;

[n,m]=size(AA);
A=AA'+AA;
A(1:n+1:end)=diag(AA);

clear Ig
clear Jg 
clear KS1 
clear BSF
clear bsparse
clear AA



% RESOLVER

% t0 = A*u';

t0_aux = load('reac_for.dat');
t0 = reshape(t0_aux(:,2:end)',3*size(node,1),1);

%%%%%%%%%%%%%%%%%%%%%%%% Loop convergence Non Linear %%%%%%%%%%%%%%%%%%%%%%

k_non = 0;

uu_trial = zeros(1,3*size(node,1)); 

TOL = 1e-2;

eps = inf;
while eps > TOL
        
uu = uu_trial;

%%%%%%%%%%%%%%%%%%%%%%%% Control Abaqus inverse %%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(strcat(filename,'_',num2str(k_non),'.inp'),'wt');
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
fprintf(fid,'*Hyperelastic, neo hooke');
fprintf(fid,'\n','');
fprintf(fid,'16.72241,6.0e-04');
fprintf(fid,'\n','');
fprintf(fid,'*STEP,NLGEOM,INC=1000');
fprintf(fid,'\n','');
fprintf(fid,'*STATIC');
fprintf(fid,'\n','');
fprintf(fid,'1.0,1.0,0.05,1.0');
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
fprintf(fid,'RF');
fprintf(fid,'\n','');
fprintf(fid,'*END STEP');
fprintf(fid,'\n','');
fprintf(fid,'*STEP,NLGEOM');
fprintf(fid,'\n','');
fprintf(fid,'*MATRIX GENERATE, STIFFNESS, ELEMENT BY ELEMENT');
fprintf(fid,'\n','');
fprintf(fid,'*MATRIX OUTPUT, STIFFNESS, FORMAT=COORDINATE');
fprintf(fid,'\n','');
fprintf(fid,'*END STEP');

fclose(fid)

command2 = strcat('abq1 & abq2 & abaqus j=',filename,'_',num2str(k_non),' user=wrt_rf.for');
status=dos(command2)


 pause(10);

% Wait for file to be created.
maxSecondsToWait = 36000; % Wait 10 hours at most.
secondsWaitedSoFar  = 0;
while secondsWaitedSoFar < maxSecondsToWait 
  if ~exist(strcat(filename,'_',num2str(k_non),'.lck'), 'file')
    break;
  end
  pause(1); % Wait 1 second.
  secondsWaitedSoFar = secondsWaitedSoFar + 1;
end


bsparse=load(strcat(filename,'_',num2str(k_non),'_STIF2.mtx'));

ind_i = bsparse(:,1)*12 + bsparse(:,2);
ind_j = bsparse(:,3);


BSF = spconvert([ind_i ind_j bsparse(:,4)]);


KS1 = (reshape(BSF',12*12,size(BSF,1)/12))';

clear Ig
clear Jg 

Ig = 3*elements(:,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0];
for i=1:11
    Ig = [Ig,3*elements(:,[1,1,1,2,2,2,3,3,3,4,4,4])-[2,1,0,2,1,0,2,1,0,2,1,0]];
end

for i=1:4
    if i==1
        Jg = 3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[2,2,2,2,2,2,2,2,2,2,2,2];
    else
        Jg = [Jg,3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[2,2,2,2,2,2,2,2,2,2,2,2]];
    end
    Jg = [Jg,3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[1,1,1,1,1,1,1,1,1,1,1,1]];
    Jg = [Jg,3*elements(:,[i,i,i,i,i,i,i,i,i,i,i,i])-[0,0,0,0,0,0,0,0,0,0,0,0]];
end


AA = sparse ( Ig , Jg , KS1, 3*size(coordinates,1) , 3*size(coordinates,1) ) ;

[n,m]=size(AA);
A=AA'+AA;
A(1:n+1:end)=diag(AA);

clear Ig
clear Jg 
clear KS1 
clear BSF
clear bsparse
clear AA

% RESOLVER

% tt_k = A*uu';

tt_aux = load('reac_for.dat');
tt_k = reshape(tt_aux(:,2:end)',3*size(node,1),1);


% %%%%%%%%%%%%%%%%%%% Problema inverso %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ind = unique([unique(elem(:,1:4));unique(fc1(:,1:3))]);


ind_a = [3*ind-2;3*ind-1;3*ind];
ind_b = setdiff(1:length(u),ind_a);

ua_s = u(ind_a)';
ub_s = u(ind_b)';

ia = length(ind_a);
ib = length(ind_b);

b_eval = -tt_k(ind_b) + A(ind_b,ind_a)*uu(ind_a)'+A(ind_b,ind_b)*uu(ind_b)';


save(strcat('b_eval_',num2str(k_non)),'b_eval');
save(strcat('tt_k_',num2str(k_non)),'tt_k');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = gmres(@(x)evaluador3(x,A,ind_a,ind_b),[ua_s; ub_s; b_eval],[],1e-3,5500);

uu_trial(ind_a) = x(1:length(ind_a));
uu_trial(ind_b) = x(length(ind_a)+1:length(ind_a)+length(ind_b));

eps = (norm(uu_trial-uu))/norm(uu_trial)

k_non = k_non+1

save(strcat('uu_trial_',num2str(k_non)),'uu_trial');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tt = tt_k;


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
fprintf(fid,'*Hyperelastic, neo hooke');
fprintf(fid,'\n','');
fprintf(fid,'16.72241,6.0e-04');
fprintf(fid,'\n','');
fprintf(fid,'*STEP,NLGEOM,INC=1000');
fprintf(fid,'\n','');
fprintf(fid,'*STATIC');
fprintf(fid,'\n','');
fprintf(fid,'1.0,1.0,0.05,1.0');
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


fprintf(fid,'*STEP,NLGEOM,INC=1000');
fprintf(fid,'\n','');
fprintf(fid,'*STATIC');
fprintf(fid,'\n','');
fprintf(fid,'1.0,1.0,0.05,1.0');
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


command2 = strcat('abq1 & abq2 & abq6141 j=',filename,'_check user=wrt_us.for');
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


