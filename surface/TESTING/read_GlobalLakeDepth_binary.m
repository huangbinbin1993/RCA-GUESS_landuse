function f=read_GlobalLakeDepth_binary

fid=fopen('/nobackup/rossby14/sm_psamu/ECOCLIMAP/RUN/GLOBAL_001/GlobalLakeDepth.001','r');

%% St Lawrence River mouth
%northp=55;
%westp=-70;
%% Sweden
northp=65;
westp=10;
%% Amazon River mouth
%northp=5;
%westp=-55;

nndeg=10;

position=(90-northp)*120*43200+(180+westp)*120;
fseek(fid,position*2,'bof');

indata=zeros(nndeg*120,nndeg*120);
for nn=1:nndeg*120
  indata(nn,:)=fread(fid,nndeg*120,'int16')';
  position=position+43200;
  fseek(fid,position*2,'bof');
end


x=westp+1/240:1/120:westp+nndeg-1/240;
y=northp-1/240:-1/120:northp-nndeg+1/240;
[XX,YY]=meshgrid(x,y);

f = indata;
if(true)
  ix=find(indata<=0);
  indata(ix)=nan;
end

%pp=pcolor(XX,YY,indata);set(pp,'edgealpha',0);colorbar('vert');
