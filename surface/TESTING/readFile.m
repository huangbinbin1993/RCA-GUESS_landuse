function [f,lon,lat]=readFile(file,marker)

fid = fopen(file,'r');
klon = fread(fid,1,'int');
klat = fread(fid,1,'int');
south = fread(fid,1,'double');
west = fread(fid,1,'double');
dlon = fread(fid,1,'double');
dlat = fread(fid,1,'double');
polon = fread(fid,1,'double');
polat = fread(fid,1,'double');

u = fread(fid,klon*klat,'double');
fclose(fid);

f=reshape(u,klon,klat);

f(find(f==marker))=nan;
%showDomain2d(polon,polat,south,west,klon,klat,dlon,dlat,file);

%visOnS2(polon,polat,south,west,klon,klat,dlon,dlat,f);


xr=west:dlon:(west+(klon-1)*dlon);
yr=south:dlat:(south+(klat-1)*dlat);
[XR,YR]=ndgrid(xr,yr);
ind = find(XR>180);
for i=ind
while XR(i)>180
   XR(i)  = XR(i)-360;
end
end

ind = find(XR<-180);
for i=ind
while XR(i)<-180
   XR(i)  = XR(i)+360;
end
end
%plot(XR,YR,'k.','MarkerSize',2)

[lon,lat]=rot2reg(XR,YR,polon,polat);

%ind = find(lon>180);
%for i=ind
%while lon(i)>180
%   lon(i)  = lon(i)-360;
%end
%end

%ind = find(lon<-180);
%for i=ind
%while lon(i)<-180
%   lon(i)  = lon(i)+360;
%end
%end


