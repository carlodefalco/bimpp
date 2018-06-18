clear all;
clc;

basename = "tumor_growth_u_";

for step = 0 : 10
    fprintf("*** Step %d ***\n", step + 1);
    t = [];
    p = [];
    u = [];
    
    for proc = 0 : 1
      fprintf("Reading and processing input file %d...\n", proc + 1);
      
      filename1 = sprintf([basename "%d_%04d"], step, proc);
      load([filename1 ".octbin.gz"]);
      t = [t, msh.t + columns(p) + 1];
      p = [p, msh.p];
      u = [u; msh.f];
    endfor
    
    msh.p = p;
    msh.t = t;
    
    fclose all;
    filename_out = sprintf([basename "_out_%d"], step);
    if (exist([filename_out ".vtu"], "file"))
        delete([filename_out ".vtu"]);
    endif
    fpl_vtk_write_field_octree(filename_out, msh,
                               {u, "u"},
                               {}, 1);
    fprintf("\n");
endfor
