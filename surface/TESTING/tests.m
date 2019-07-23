%clear;close all

f=read_GlobalLakeDepth_binary;

[n,m] = size(f);

N = 40;%this is how many points are to be used in the integral

D = zeros(floor(n/N),floor(m/N));
M = zeros(floor(n/N),floor(m/N));
S = zeros(floor(n/N),floor(m/N));
[nn,mm] = size(D);
shallow = 0;
medium = 3;
deep = 8;%12;
for j=1:mm
    for i=1:nn
        g = f(((i-1)*N +1):((i-1)*N +N) ,((j-1)*N +1):((j-1)*N+N))/10;
        %from decimeter to meters
        nnz = find(g>0);
        
        if(length(nnz)>0)
            h = g(nnz);
            
            d = find(h>deep);
            if(length(d)>0)
                dd = h(d);
                D(i,j) = sum(dd(:))/(length(d));
            else
                D(i,j) = nan;
            end
            
            med = find(h>medium & h<=deep);
            if(length(med)>0)
                medmed = h(med);
                M(i,j) = sum(medmed(:))/(length(med));
            else
               M(i,j) = nan;
            end
            
            s = find(h>shallow & h<=medium);
            if(length(s)>0)
                ss = h(s);
                S(i,j) = sum(ss(:))/(length(s));
            else
                S(i,j) = nan;
            end
        else
            D(i,j) = nan;
            M(i,j) = nan;
            S(i,j) = nan;
        end
    end
end

%close all

figure
S(find(S==0))=nan;
contourf(flipud(S))
title('shallow')
colorbar
axis equal
figure
M(find(M==0))=nan;
contourf(flipud(M))
title('medium')
colorbar
axis equal
figure
D(find(D==0))=nan;
contourf(flipud(D))
title('deep')
colorbar
axis equal

% $$$ figure
% $$$ f(find(f==0))=nan;
% $$$ contourf(flipud(f/10))
% $$$ colorbar
% $$$ axis equal                              % 