#!/usr/bin/env octave 
pkg load bim

xx = linspace (0,8,120);
yy = linspace (0,1,11);
zz = linspace (0,1,15);
msh = msh3m_structured_mesh (xx, yy, zz, 1, 1:6);

msh.t(1:4, :) -= 1;
msh.e(1:3, :) -= 1;

fid = fopen ("mesh_in.msh", 'w');

fprintf (fid, "%d %d %d \n\n", columns (msh.p), columns (msh.t), columns (msh.e));

fprintf (fid, "%17.17g ", msh.p(:));
fprintf (fid, "\n\n");

fprintf (fid, "%d ", msh.t(:));
fprintf (fid, "\n\n");

fprintf (fid, "%d ", msh.e(:));
fprintf (fid, "\n");

fclose (fid);

system ("make");
system ("./bim3_test");

mesh_out

msh2.p =p;
msh2.t =t;
msh2.e =e;

delete ("mesh_out.vtu")
fpl_vtk_write_field ("mesh_out", msh2, {p(1,:)', 'p1'}, {}, 1);
