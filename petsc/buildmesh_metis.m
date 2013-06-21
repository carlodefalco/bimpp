function buildmesh_metis(p,e,t,data,nproc)

% buildmesh_metis(p,e,t,data,nproc)


nelem = size(t,2);

f=fopen('t','w');
fprintf(f,'%d %d\n',nelem,1);
fprintf(f,'%d %d %d\n',t(1:3,:));
fclose(f);

f=fopen('p','w');
fprintf(f,'%d %d \n',p);
fclose(f);

f=fopen('data','w');
fprintf(f,'%d %d %d\n',data);
fclose(f);

bcnodes=sort(unique([e(1,:) e(2,:)]));
f=fopen('bc','w');
fprintf(f,'%d \n',bcnodes);
fclose(f);

system(['$METIS/partnmesh t ' num2str(nproc)])
