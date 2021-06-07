[x_read2,y_read2,z_read2] = ndgrid((1:size(dispField.X,1))*Lx,(1:size(dispField.X,2))*Ly,(1:size(dispField.X,3))*Lz);




tmpX = griddedInterpolant(x_read2,y_read2,z_read2,dispField.Y);
tmpY = griddedInterpolant(x_read2,y_read2,z_read2,dispField.X);

tmp_Y = tmpY(node(:,1),node(:,2),node(:,3));
tmp_X = tmpX(node(:,1),node(:,2),node(:,3));

tic
tmpZ = griddedInterpolant(x_read2,y_read2,z_read2,dispField.Z);
tmp_Z = tmpZ(node(:,1),node(:,2),node(:,3));
disp(['Time interpolating using gridded interpolant: ', num2str(toc)])

tic
FUZ = scatteredInterpolant(x_read',y_read',z_read',double(w_read));
UZ = FUZ(node(:,1),node(:,2),node(:,3));
disp(['Time interpolating using scattered interpolant: ', num2str(toc)]);