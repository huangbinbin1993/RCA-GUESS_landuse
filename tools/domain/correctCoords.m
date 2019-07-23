function xOut = correctCoords(x)
xOut = x;
[m,n]=size(x);

xOut(1) = xOut(1)+pi;
xOut(2) = xOut(2)+pi/2;

xOut(1) = mod(xOut(1),2*pi);
xOut(2) = mod(xOut(2),2*pi);
if(xOut(2) > pi)
    xOut(2) = pi - xOut(2);
    xOut(1) = xOut(1)+pi;
end


 xOut(1) = mod(xOut(1),2*pi);   
    
%if(abs(xOut(2)) > pi)
%    xOut(2) = mod(xOut(2),pi); %change its interval to 0,2pi
%    xOut(1) = xOut(1)+pi;
%end

%if(abs(xOut(1))>2*pi)
%    xOut(1) = mod(xOut(1),2*pi);
%end

xOut(1) = xOut(1)-pi;
xOut(2) = xOut(2)-pi/2;

