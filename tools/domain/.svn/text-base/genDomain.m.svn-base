function [lon,lat]=genDomain(polon,polat,south,west,klon,klat,dlon,dlat)
%genDomain(polon,polat,south,west,klon,klat,dlon,dlat)
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
