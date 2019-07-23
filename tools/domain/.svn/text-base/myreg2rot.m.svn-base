function [xrot,J] = myreg2rot(xreg, polreg);


xrot = 100*ones(size(xreg));
J=zeros(2,2);
%sin(x+pi/2)=cos(x)
%cos(x+pi/2)=-sin(x)

polrot(1) = 0;
polrot(2) = -pi/2;

polreg = correctCoords(polreg);
polon = polreg(1);
polat = polreg(2);

xrot(2) = asin( -sin(polat).*sin(xreg(2)) - cos(polat).*cos(xreg(2)).*cos(xreg(1)-polon));

J(2,1)=cos(xreg(2))*cos(polat)*sin(polon - xreg(1))/sqrt(1 - ...
                                                  (cos(xreg(2))*cos(polat)*cos(polon - xreg(1))...
                                                  + sin(xreg(2))*sin(polat))^ 2 );

J(2,2)=(cos(xreg(2))*sin(polat)*cos(polon - xreg(1)) -sin(xreg(2))*cos(polat))/...
       sqrt(1 - (cos(xreg(2))*cos(polat)*cos(polon - xreg(1)) + sin(xreg(2))*sin(polat))^ 2 );


sxrot= cos(xreg(2)).*sin(xreg(1)-polon)./cos(xrot(2));

xrot(1)= sign(sxrot)*...
           acos((-sin(polat).*cos(xreg(2)).*cos((xreg(1)-polon)) + ...
                 cos(polat).*sin(xreg(2)))./cos(xrot(2)));

J(1,1)=sign(sxrot)*(sin(polat)*sin(xreg(1) - polon)*cos(xreg(2)))/...
       (sqrt(1- ((cos(polat)*sin(xreg(2)) - sin(polat)*cos(xreg(1) - polon)*cos(xreg(2)))^2)/...
             (cos(xrot(2)))^2 ));

J(1,2) = - (sign(sxrot)*(- sin(polat)*sin(xreg(2)) - cos(polat)*cos(xreg(1) - polon)*cos(xreg(2))))/...
         (sqrt(1 - ((cos(polat)*sin(xreg(2)) - sin(polat)*cos(xreg(1) - polon)*cos(xreg(2)))^2)/...
               (cos(xrot(2)))^2));

xrot = correctCoords(xrot);