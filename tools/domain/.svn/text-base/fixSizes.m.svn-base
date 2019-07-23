function [klon,klat] = fixSizes(klon,klat)

if mod(klon,2) ~=0
  klon = klon+1;
end
if mod(klat,2) ~=0
  klat = klat+1;
end

N = [klon klat] -6;

for i=1:length(N)
  primes = factor(N(i));
  primes(find(primes==2))=[];
  primes(find(primes==3))=[];
  primes(find(primes==5))=[];
  while ~isempty(primes)
    N(i)=N(i)+1;
    primes = factor(N(i));
    primes(find(primes==2))=[];
    primes(find(primes==3))=[];
    primes(find(primes==5))=[];
  end
  
end

N = N+6;
klon = N(1);
klat = N(2);