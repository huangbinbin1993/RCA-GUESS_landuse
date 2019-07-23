function [f,x,y]=readTESTFile(file)
file
fid = fopen(file,'r')
klon = fread(fid,1,'int')
klat = fread(fid,1,'int')
south = fread(fid,1,'double')
west = fread(fid,1,'double')
dlon = fread(fid,1,'double')
dlat = fread(fid,1,'double')
polon = fread(fid,1,'double')
polat = fread(fid,1,'double')

u = fread(fid,klon*klat,'double');
f=reshape(u,klon,klat);

  u = fread(fid,klon*klat,'double');
  x=reshape(u,klon,klat);
  u = fread(fid,klon*klat,'double');
  y=reshape(u,klon,klat);

fclose(fid);

