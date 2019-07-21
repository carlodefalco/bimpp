function [] = Export_ch_sol(n_processors = 1,range = 10000)
% Export the solution


% Go from the current directory to script
addpath (canonicalize_file_name ("../../script/m/"));


export_tmesh_data ("cahn_hilliard_u_00_%4.4d_%4.4d", {"cahn_hilliard_u_00_%4.4d_%4.4d", "cahn_hilliard_u_01_%4.4d_%4.4d"}, {"u_00","u_01"},{},{},"solution",0:range,n_processors);

%export_tmesh_data ("cahn_hilliard_u_00_%4.4d_%4.4d", {"cahn_hilliard_u_00_%4.4d_%4.4d", "cahn_hilliard_u_01_%4.4d_%4.4d", "cahn_hilliard_u_02_%4.4d_%4.4d","cahn_hilliard_u_03_%4.4d_%4.4d"}, {"u_00","u_01", "u_02", "u_03"},{},{},"solution",0:range,n_processors);

end
