%REG2ROT converts regular coordinates to rotated ones.
%
%   [X_ROT, Y_ROT]= REG2ROT(X_REG, Y_REG, X_C, Y_C)
%
%   returns the rotated coordinates X_ROT and Y_ROT corresponding to
%   the regular coordinates X_REG and Y_REG. The regular coordinates
%   of the southpole in the rotated grid is given by X_C, Y_C.

function [xrot, yrot]= reg2rot(xreg, yreg, xc, yc);

rad= pi/180;
syc= sin(rad * (yc+90));
cyc= cos(rad * (yc+90));

xmxc= rad * (xreg - xc);
sxmxc= sin(xmxc);
cxmxc= cos(xmxc);
syreg= sin(rad * yreg);
cyreg= cos(rad * yreg);
syrot= cyc.*syreg - syc.*cyreg.*cxmxc;
ix=find(isnan(syrot)==1);
syrot= min(max(syrot, -1), 1);
syrot(ix)=ones(length(ix),1).*nan;

yrot= 1/rad * asin(syrot);

cyrot= cos(rad * yrot);
cxrot= (cyc.*cyreg.*cxmxc + syc.*syreg) ./ cyrot;
ix=find(isnan(cxrot)==1);
cxrot= min(max(cxrot, -1), 1);
cxrot(ix)=ones(length(ix),1).*nan;
sxrot= cyreg.*sxmxc./cyrot;

xrot= 1/rad * acos(cxrot);

tmp= sxrot < 0;
xrot= xrot.*(1 - 2.*tmp);
