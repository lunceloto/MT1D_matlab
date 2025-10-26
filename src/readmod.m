function [dx,dy,dz,rho]=readmod(modname)
fid=fopen(modname);
if fid==-1
    disp("error:Can't open mod");
end
fgetl(fid);
line = fgetl(fid);
[grids] = sscanf(line,'%d',[4 1]);
rho=zeros(grids(1),grids(2),grids(3));
dx=fscanf(fid,'%f',grids(1));
dy=fscanf(fid,'%f',grids(2));
dz=fscanf(fid,'%f',grids(3));
fgetl(fid);
for i=1:grids(3)
    fgetl(fid);
    for j=1:grids(2)
        rho(:,j,i)=fscanf(fid,'%f',grids(2));
    end
end
fclose(fid);
end

