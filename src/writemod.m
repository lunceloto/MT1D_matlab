function writemod(dx,dy,dz,rho,fname)
fid = fopen(fname,'w+');
char=strcat('#',fname);
fprintf(fid,'%s\n',char);
grids=size(rho);
fprintf(fid,'%5d %5d %5d %5d',[grids(1),grids(2),grids(3),0]);
fprintf(fid,'%7s\n',' LINEAR');
fprintf(fid,'%13.5f\t',dx);fprintf(fid,'\n');
fprintf(fid,'%13.5f\t',dy);fprintf(fid,'\n');
fprintf(fid,'%13.5f\t',dz);fprintf(fid,'\n');
for i= 1:grids(3)
    fprintf(fid,'\n');
    for j = 1:grids(2)
        fprintf(fid,'%13.5f\t',rho(:,j,i));fprintf(fid,'\n');
    end
end
fclose(fid);
end
