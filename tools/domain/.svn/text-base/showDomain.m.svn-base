function showDomain(polon,polat,south,west,klon,klat,dlon,dlat)
%showDomain(polon,polat,south,west,klon,klat,dlon,dlat)
% GRID
east = west + (klon-1)*dlon;
north = south +(klat-1)*dlat;
xr=linspace(west,east,klon); %deg
yr=linspace(south,north,klat); %deg

[XR,YR]=ndgrid(xr,yr); %deg
[x,y]=rot2reg(XR,YR,polon,polat); %deg
x=x*pi/180;
y=y*pi/180;
x=x+pi;
y=y+pi/2;

rr=1.01;

Xgrid = rr.*sin(y).*sin(x);
Ygrid = rr*sin(y).*cos(x);
Zgrid = rr.*cos(y);
%end grid

%doing coastal lines
load n_coast;

n_coast = n_coast*pi/180;
lon = n_coast(:,1)+pi;
lat = n_coast(:,2)+pi/2;

Xmap = rr.*sin(lat).*sin(lon);
Ymap = rr.*sin(lat).*cos(lon);
Zmap = rr.*cos(lat);
%end coast

%figure(1);
hold on
plot3(Xmap,Ymap,Zmap,'k')
plot3(Xgrid,Ygrid,Zgrid,'k.','MarkerSize',4)


%a blue globe
theta = linspace(0,2*pi,200);%lon
phi = linspace(0,pi,100); %lat
r = linspace(0.9,1.1);
[T,P,R]=ndgrid(theta,phi,r);

X = R.*sin(P).*cos(T);
Y = R.*sin(P).*sin(T);
Z = R.*cos(P);
globe = X.^2 + Y.^2 + Z.^2;

p=patch(isosurface(X,Y,Z,globe,1.0));
set(p, 'FaceColor', 'blue', 'EdgeColor', 'none');


%vis options
daspect([1 1 1]); axis tight; 
camlight; lighting none;
hold on
axis vis3d
view(polon-180,271-polat)