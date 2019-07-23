function [polon,polat,south,west,klon,klat,dlon,dlat,name]=Cordex(region)

dlon=0.44;
dlat=0.44;
switch region
  case 1
    name ='South_America';
    RotPole=[-56.06;70.6]; 
    TLC= [144.14;34.54];
    Nx=146;
    Ny=167;
  case 2
    name='Central_America';
    RotPole=[113.98;75.74];
    TLC=[307.20;20.68];
    Nx=210;
    Ny=113;
  case 3
    name='North_America';
    RotPole=[83.0;42.5];
    TLC=[326.12;28.36];
    Nx=155;
    Ny=130;
  case 4
    name='Europe';
    RotPole=[198.0;39.25];
    TLC=[331.79;21.67];
    Nx=106;
    Ny=103;

  case 5
    name='Africa';
    RotPole=[180.0;90.0];
    TLC=[335.36;42.24];
    Nx=194;
    Ny=201;
    
  case 6
    name='West_Asia';
    RotPole=[236.66;79.95];
    TLC=[327.88;35.20];
    Nx=193;
    Ny=130;
    
  case 7
    name='East_Asia';
    RotPole=[295.22;77.61];
    TLC=[319.08;46.20];
    Nx=203;
    Ny=167;
    
  case 8
    name='Central_Asia';
    RotPole=[256.61;43.48];
    TLC=[325.68;22.88];
    Nx=153;
    Ny=100;
    
  case 9
    name='Australasia';
    RotPole=[321.28;-60.31];
    TLC=[142.16;33.44];
    Nx=200;
    Ny=129;
    
  case 10
    name='Antarctica';
    RotPole=[-166.92;6.08]; %[13.08;-6.08];
    TLC=[152.94;14.3];
    Nx=125;
    Ny=97;
    
  case 11
    name='Arctic';
    RotPole=[0.0;6.55];
    TLC=[337.12;33.88];
    Nx=116;
    Ny=133;
    
  case 12
    name='Mediterranean_domain';
    RotPole=[198.0;39.25];
    TLC=[-23.0;5.72];%[336.78;5.94];
    Nx=98;
    Ny=63;
    
  case 13
    name='High_Res_Europe';
    RotPole=[198.0;39.25];
    TLC=[331.79;21.67];
    Nx=4*106;
    Ny=4*103;
    dlon=0.11;
    dlat=0.11;
    
    case 14
    name='High_Res_Arctic';
    RotPole=[0.0;6.55];
    TLC=[337.12;33.88];
    Nx=2*116;
    Ny=2*133;
    dlon=0.22;
    dlat=0.22;
  case 15
    name='High_Res_Africa';
    RotPole=[180.0;90.0];
    TLC=[335.36;42.24];
    Nx=2*194;
    Ny=2*201;
    dlon=0.22;
    dlat=0.22;
  otherwise
    name='Unknown Cordex Domain';
    
end

%do conversion to RCA domain
polon = RotPole(1); 
polat =  RotPole(2)+180.0; %Cordex specifies the North Pole
if(abs(polat)>90)
    if(polat>90)
        polat = 90 - (polat-90);
        polon = polon +180;
    else
        polat = -90 - (polat+90);
        polon = polon +180;
    end
end
while(abs(polon)>180)
    polon = polon - sign(polon)*360;
end

n_relax_points = 8;
npassive = 2;

klon = Nx;
klat = Ny;
BLC = [TLC(1);TLC(2)-(Ny-1)*dlat];

west = BLC(1)-(n_relax_points+npassive)*dlon; 
south = BLC(2)-(n_relax_points+npassive)*dlat;
if(abs(south)>90)
    if(south>90)
        south = 90 - (south-90);
        west = west +180;
    else
        south = -90 - (south+90);
        west = west +180;
    end
end
while(abs(west)>180)
    west = west - sign(west)*360;
end
klon = klon + 2*(n_relax_points+npassive);
klat = klat + 2*(n_relax_points+npassive);


klon = fixRCAFFTsize(klon);
klat = fixRCAFFTsize(klat);



disp('&domain');
disp(['klon_global = ',int2str(klon)]);
disp(['klat_global = ',int2str(klat)]);
disp(['klev_global = ',int2str(40)]);
disp(['south = ',num2str(south)]);
disp(['west = ',num2str(west)]);
disp(['dlon = ',num2str(dlon)]);
disp(['dlat = ',num2str(dlat)]);
disp(['polon =',num2str(polon)]);
disp(['polat = ',num2str(polat)]);



disp(['&nampos']);
disp(['nppstr = 5']);
disp(['lunlis= 57']);
disp(['lundip= 71']);
disp(['lprint=.false.']);
disp(['iprint= 10']);
disp(['jprint= 10']);
disp(['lomega=.true.']);
disp(['hisnam=''fc''']);
disp(['npplon=',int2str(Nx)]);
disp(['npplat=',int2str(Ny)]);
disp(['iminpp=',int2str(n_relax_points+npassive+1)]);
disp(['jminpp=',int2str(n_relax_points+npassive+1)]);




% $$$ if(false)
% $$$     path=['/home/sm_marku/rca/reference_domains/Cordex/' name];
% $$$     name = [path '/namelists.dat'];
% $$$     fid = fopen(name,'w');
% $$$     fprintf(fid,'&domain\n');
% $$$     fprintf(fid,'klon_global = %d\n',klon);
% $$$     fprintf(fid,'klat_global = %d\n',klat);
% $$$     fprintf(fid,'klev_global = %d\n',40);
% $$$     fprintf(fid,'south = %12.8f\n',south);
% $$$     fprintf(fid,'west = %12.8f\n',west);
% $$$     fprintf(fid,'dlon = %12.8f\n',dlon);
% $$$     fprintf(fid,'dlat = %12.8f\n',dlat);
% $$$     fprintf(fid,'polon = %12.8f\n',polon);
% $$$     fprintf(fid,'polat = %12.8f\n',polat);
% $$$     fprintf(fid,'/\n\n');
% $$$     
% $$$     
% $$$     fprintf(fid,'&nampos\n');
% $$$     fprintf(fid,'nppstr = 5\n');
% $$$     fprintf(fid,'lunlis= 57\n');
% $$$     fprintf(fid,'lundip= 71\n');
% $$$     fprintf(fid,'lposton=.false. \n');
% $$$     fprintf(fid,'lprint=.false.\n');
% $$$     fprintf(fid,'iprint= 10\n');
% $$$     fprintf(fid,'jprint= 10\n');
% $$$     fprintf(fid,'lomega=.true.\n');
% $$$     fprintf(fid,'hisnam=''fc''\n');
% $$$     fprintf(fid,'ntimep=1\n');
% $$$     fprintf(fid,'npplon=%d\n',Nx);
% $$$     fprintf(fid,'npplat=%d\n',Ny);
% $$$     fprintf(fid,'iminpp=%d\n',n_relax_points+npassive+1);
% $$$     fprintf(fid,'jminpp=%d\n',n_relax_points+npassive+1);
% $$$     fprintf(fid,'/\n\n');
% $$$     fclose(fid);
% $$$ end
