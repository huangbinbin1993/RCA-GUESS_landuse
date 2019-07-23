lonN = [-180:40:180];
lonS = [-180:60:180];
lat = [-90 -60 -10 40 90];
hold on

for k=1:length(lat)
    if(k<=2)
        plot(lonS,lat(k)*ones(size(lonS)))
    else
        plot(lonN,lat(k)*ones(size(lonN)))
    end
end

for k=1:length(lonN)
    plot(lonN(k)*ones(size(lat(2:end))),lat(2:end))
end
for k=1:length(lonS)
    plot(lonS(k)*ones(size(lat(1:2))),lat(1:2))
end
    
