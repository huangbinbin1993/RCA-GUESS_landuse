


for k=[1 2 3 4 5 6 7 8 9 12]
%for k=10
        [polon,polat,south,west,klon,klat,dlon,dlat,name]=Cordex(k);
        showDomain2d(polon,polat,south,west,klon,klat,dlon,dlat,name);
        hold on;
end
axis([-180 180 -90 90])

plotGtopo();