clear;close all;
disp('This program will help you to choose the domain parameters');
dlon = input('Input dlon:  ');
dlat = input('Input dlat:  ');
key = input('Do you wish to input graphically? yes=1, no=0:   ');
%dlon = 0.25;dlat=0.25;key=0;
%iterativeHH = input(['Are you going to use the interative HH-solver ' ...
%                    '(not supported currently) yes=1, no=0::']);
figure(1);
load n_coast;
hold on
plot(n_coast(:,1) ,n_coast(:,2) ) %deg

xI=[];
yI=[];
if(key==1)
    disp(['Click four corners on the map. This will be the bounding box ' ...
      'of your computational domain'])
    for i=1:4
       [xt,yt] = ginput(1);
       xI=[xI;xt]; %deg
       yI=[yI;yt]; %deg
       plot(xI,yI,'ko-')
    end
else
  xinmin = input('Input lon_min: ');
  yinmin = input('Input lat_min: ');
  plot(xinmin,yinmin,'o')
  xinmax = input('Input lon_max: ');
  yinmax = input('Input lat_max: ');
  plot(xinmax,yinmax,'o')
  xI = [xinmin xinmin xinmax xinmax]';
  yI = [yinmin yinmax yinmax yinmin]';
end

xI=xI;
yI=yI;
xp = [xI; xI(1)];%deg
yp = [yI; yI(1)]; %deg
plot(xp,yp,'ko-','MarkerSize',10); %bbox border

xm = [mean(xI); mean(yI)]; %deg
xmin = min(xI); %deg
ymin = min(yI); %deg
xmax = max(xI); %deg
ymax = max(yI); %deg
plot(xm(1),xm(2),'kp-','MarkerSize',10)%midpoint
drawnow

xtmp = [xmin xmin xmax xmax xmin];
ytmp = [ymin ymax ymax ymin ymin];

pol = [0.5; 0];;%starting guess in Rad
              %[pol(1),pol(2)] = ginput(1);
              %pol = xm+100*eps;
              %pol = pol*pi/180;

plot(0,0,'rp-','MarkerSize',10)
x = [1; 1];%starting guess in Rad
iter = 0;
maxIter = 10000;
N=zeros(maxIter+1);
disp('  iter   norm   polon  polat  ')
while (norm(x,inf)>0.001 & iter<maxIter)
    [x,J] = myreg2rot(xm*pi/180,pol); %rad
    x = correctCoords(x); %rad
    pol = pol - J\x; %rad
%    if(norm(x,inf)<0.003)
%        disp([iter norm(x,inf) pol'*180/pi])
%    end
    
    N(iter+1)=norm(x);
    iter = iter+1;
end
%disp(' iter        norm')
%disp([iter norm(x,inf)])

pol = correctCoords(pol);
pol = pol*180/pi; %deg


[XX,YY] = reg2rot(xp,yp,pol(1),pol(2));
west = min(XX);east = max(XX);south = min(YY);north = max(YY);
%plot([west west east east],[south north north south],'r*') %deg


%Dlon = max([abs(east) abs(west)]); Dlat = max([abs(north)
%abs(south)]); %deg
Dlon = east-west; Dlat = north-south; %deg
%east=Dlon;west=-Dlon;north=Dlat;south=-Dlat;
klon = floor(Dlon/dlon)+1; %deg
klat = floor(Dlat/dlat)+1; %deg

%check if klon,klat are valid for the fft-based helmholtz-solver
%and shift accordingly

if(true) %when we are using the fft-based helmholtz-solver there
         %are conditions on the sizes
  [klon,klat] = fixSizes(klon,klat);
  north = south + (klat-1)*dlat;
  east = west + (klon-1)*dlon;
end

%plot([west west east east],[south north north south],'r*') %deg
%plot([xmin xmin xmax xmax],[ymin ymax ymax ymin],'k*')  %deg

xr=linspace(west,east,klon); %deg
yr=linspace(south,north,klat); %deg

[XR,YR]=ndgrid(xr,yr); %deg
[n,m]=size(XR);
xx=XR;yy=YR;
[x,y]=rot2reg(xx,yy,pol(1),pol(2)); %deg

for i=1:length(x(:))
    tmp = correctCoords([x(i)*pi/180;y(i)*pi/180]); %red
    x(i) = tmp(1)*180/pi;%deg
    y(i) = tmp(2)*180/pi; %deg
end
for i=1:length(xx(:))
    tmp = correctCoords([xx(i)*pi/180;yy(i)*pi/180]); %red
    xx(i) = tmp(1)*180/pi;%deg
    yx(i) = tmp(2)*180/pi; %deg
end

plot(x,y,'k.','MarkerSize',2)
%plot(xx,yy,'r.','MarkerSize',2);


polon = pol(1);
polat = pol(2);

disp('   polon     polat     south     west      klon_global klat_global')
disp([polon polat south west klon klat])
