function export_tmesh_data (msh_basename,
                            nodedata_basenames, nodedata_names,
                            celldata_basenames, celldata_names,
                            out_name = "out",
                            steps = 0:10, nprocs = 1)

  nnodefields = numel (nodedata_names);
  ncellfields = numel (celldata_names);
  
  if (nnodefields == 0)
    n = {};
    nodedata_names = {};
  endif
  
  if (ncellfields == 0)
    c = {};
    celldata_names = {};
  endif
  
  for step = steps
    fprintf("*** Step %d ***\n", step);
    t = [];
    p = [];

    for ii = 1 : nnodefields
      n{ii} = [];
    endfor

    for ii = 1 : ncellfields
      c{ii} = [];
    endfor
    
    for proc = 0 : nprocs - 1
      fprintf("Reading and processing input file %d...\n", proc);
      
      filename = sprintf (msh_basename, step, proc);
      load ([filename ".octbin.gz"]);
      
      t = [t, msh.t + columns(p) + 1];
      p = [p, msh.p];

      for ii = 1 : nnodefields
        filename = sprintf (nodedata_basenames{ii}, step, proc);
        load ([filename ".octbin.gz"]);
        n{ii} = [n{ii}; msh.f];
      endfor

      for ii = 1 : ncellfields
        filename = sprintf (celldata_basenames{ii}, step, proc);
        load ([filename ".octbin.gz"]);
        c{ii} = [c{ii}; msh.f];
      endfor
    endfor
    
    msh.p = p;
    msh.t = t;
    
    filename_out = sprintf ([out_name "_%d"], step);
    if (exist ([filename_out ".vtu"], "file"))
      delete ([filename_out ".vtu"]);
    endif

    fpl_vtk_write_field_octree (filename_out, msh,
                                [n; nodedata_names].',
                                [c; celldata_names].', 1);
    
    fprintf("\n");
  endfor
endfunction
