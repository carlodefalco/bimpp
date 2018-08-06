function export_tmesh_data (msh_basename, data_basenames,
                            data_names, out_name = "out",
                            steps = 0:10, nprocs = 1)

  nfields = numel (data_names);
  for step = steps
    
    fprintf("*** Step %d ***\n", step + 1);
    t = [];
    p = [];

    for ii = 1 : numel (nfields)
      u{ii} = [];
    endfor
    
    for proc = 0 : nprocs
      
      fprintf("Reading and processing input file %d...\n", proc + 1);
      
      filename1 = sprintf([msh_basename "%d_%04d"], step, proc);
      load ([filename1 ".octbin.gz"]);
      
      t = [t, msh.t + columns(p) + 1];
      p = [p, msh.p];

      for ii = 1 : numel (nfields)
        filename1 = sprintf([data_basenames{ii} "%d_%04d"], step, proc);
        load ([filename1 ".octbin.gz"]);
        u{ii} = [u{ii}; msh.f];
      endfor
      
    endfor
    
    msh.p = p;
    msh.t = t;
    
    filename_out = sprintf ([basename "_out_%d"], step);
    if (exist([filename_out ".vtu"], "file"))
      delete([filename_out ".vtu"]);
    endif
    
    fpl_vtk_write_field_octree (filename_out, msh,
                                [u; data_names{:}].',
                                {}, 1);
    
    fprintf("\n");

  endfor

endfunction
