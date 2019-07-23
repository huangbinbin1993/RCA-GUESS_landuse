%ROT2REG converts rotated coordinates to regular ones.
%
%   [X_REG, Y_REG]= ROT2REG(X_ROT, Y_ROT, X_C, Y_C)
%
%   returns the regular coordinates X_REG and Y_REG corresponding to
%   the rotated coordinates X_ROT and Y_ROT. The regular coordinates
%   of the southpole in the rotated grid is given by X_C, Y_C.

function [xreg, yreg]= rot2reg(xrot, yrot, xc, yc);

rad= pi/180;
syc= sin(rad * (yc+90));
cyc= cos(rad * (yc+90));

sxrot= sin(rad * xrot);
cxrot= cos(rad * xrot);
syrot= sin(rad * yrot);
cyrot= cos(rad * yrot);
syreg= cyc.*syrot + syc.*cyrot.*cxrot;
syreg= min(max(syreg, -1), 1);

yreg= 1/rad * asin(syreg);

cyreg= cos(rad * yreg);
cxmxc= (cyc.*cyrot.*cxrot - syc.*syrot)./cyreg;
cxmxc= min(max(cxmxc, -1), 1);
sxmxc= cyrot.*sxrot./cyreg;
xmxc= 1/rad * acos(cxmxc);

tmp= sxmxc < 0;
xmxc= xmxc.*(1 - 2*tmp);;
xreg= xmxc + xc;