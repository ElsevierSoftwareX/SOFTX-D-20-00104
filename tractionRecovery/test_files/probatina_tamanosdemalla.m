

load([data_path filesep file_name],'cellMask');
cellMask.mask3D = cellMask.mask3D(1:300,:,:);


maxvols_gel = 100:100:600;
maxsizes_cell = 1:0.5:5;

container = zeros(length(maxvols_gel),2*length(maxsizes_cell));

for ii = 1:length(maxvols_gel)
    disp(['maxvol ' num2str(ii) '/' num2str(length(maxvols_gel))]);
    c=1;
    for jj=1:2:2*length(maxsizes_cell)
        disp(['...maxsize ' num2str(jj) '/' num2str(length(maxsizes_cell))]);
        maxvol_gel = maxvols_gel(ii); % How fine we mesh the gel. Something between 200 and 500
        maxsize_cell = maxsizes_cell(c); % How fine we mesh the sprout surface
        c=c+1;
        [~,elem_gel,elem_cell,~] = getNodesMask(cellMask.mask3D,Lx,Ly,Lz,maxsize_cell,maxvol_gel,write_options);
        
        container(ii,jj:jj+1) = [size(elem_gel,1),size(elem_cell,1)];
    end
end
disp('DONE!');


figure
subplot(131)
surf(maxsizes_cell,maxvols_gel,container(:,1:2:end));
xlabel('maxsizes cell');ylabel('maxvols gel'); title('Elements gel');
view([0 90]);
colormap(jet(256));colorbar;
subplot(132)
surf(maxsizes_cell,maxvols_gel,container(:,2:2:end));
xlabel('maxsizes cell');title('Elements cell');
view([0 90]);
colormap(jet(256));colorbar;
total_elements = container(:,1:2:end)+container(:,2:2:end);
subplot(133)
surf(maxsizes_cell,maxvols_gel,total_elements);
xlabel('maxsizes cell');title('TOTAL');
view([0 90]);
colormap(jet(256));colorbar;

