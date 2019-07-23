function showDomain2d(polon,polat,south,west,klon,klat,dlon,dlat,name)
%showDomain2d(polon,polat,south,west,klon,klat,dlon,dlat)
% GRID
east = west + (klon-1)*dlon;
north = south +(klat-1)*dlat;
xr=linspace(west,east,klon); %deg
yr=linspace(south,north,klat); %deg

[XR,YR]=ndgrid(xr,yr); %deg
[lon,lat]=rot2reg(XR,YR,polon,polat); %deg

ind =find(lat>90);
lat(ind) = 90 - (lat(ind)-90);
lon(ind) = lon(ind) +180;

ind =find(lat<-90);
lat(ind) = -90 - (lat(ind)+90);
lon(ind) = lon(ind) +180;

ind = find(abs(lon)>180);
while(~isempty(ind))
    lon(ind) = lon(ind) - sign(lon(ind))*360;
    ind = find(abs(lon)>180);
end

%doing coastal lines
load n_coast;
%keyboard
xb = [lon(1:klon,1); (lon(klon,2:klat))'; flipud(lon(1:klon-1,klat)); (fliplr(lon(1,1:klat-1)))'];
yb = [lat(1:klon,1); (lat(klon,2:klat))'; flipud(lat(1:klon-1,klat)); (fliplr(lat(1,1:klat-1)))'];

plot(n_coast(:,1),n_coast(:,2),xb,yb,'*')%lon,lat,'k.','markerSize',4)
axis([min(lon(:)) max(lon(:)) min(lat(:)) max(lat(:))]) 
title(name)
xlabel('LON');ylabel('LAT')
print('-depsc',([name '.eps']));