function f=fixRCAFFTsize(N)
    
    N = N-6;
    if(mod(N,2)~=0)
        N = N+1; %make even!
    end
    prims = factor(N);
    prims(find(prims==2 |prims==3 | prims==5))=[];
    while ~isempty(prims)
        N = N+2;
        prims = factor(N);
        prims(find(prims==2 |prims==3 | prims==5))=[];
    end
    
    f = N+6;
    