function vizBox(s,n,w,e)
N = 3;

x = [s s n n s];
y = [w e e w w];

load n_coast;


plot(n_coast(:,1),n_coast(:,2),x,y)