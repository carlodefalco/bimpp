addpath(pwd);

%[p,e,t]=initmesh('squareg','Hmax',1/5);
[p,e,t]=poimesh('squareg',100,100);
p = (p+1)/2;

x=p(1,:);
y=p(2,:);

u = x.*0;
cerchio = find(sqrt((x-.3).^2+(y-.7).^2)<.1);
u(cerchio) = 1;

f = x.*0;
v = 10*(x-y);

data = [v;f;u];

for nproc=[2,4,8,16,32]
  system(['mkdir nproc' num2str(nproc)]);
  system(['cp plotseries.m nproc' num2str(nproc)]);
  cd(['nproc' num2str(nproc)]);
  save matdata p e t data
  buildmesh_metis(p,e,t,data,nproc);
  cd ..
end
