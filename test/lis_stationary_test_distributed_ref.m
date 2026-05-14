% Reference solution for lis_stationary_test_distributed.cpp
% Assembles the same 55-node 2D quad mesh and solves with Octave's backslash.

N   = 55;  % nodes: 11 columns x 5 rows
nel = 40;  % elements: 10 columns x 4 rows

A = sparse (N, N);
b = ones (N, 1);

locmatrix = [ 2, -1, -1,  0;
             -1,  2,  0, -1;
             -1,  0,  2, -1;
              0, -1, -1,  2];

% Mirrors the C++ conn(inode, iel) function (0-based in, +1 for Octave)
function g = conn (inode, iel)
  elcol = floor (iel / 4);
  elrow = mod (iel, 4);
  g = elcol * 5 + elrow;
  if inode == 1, g += 1; end
  if inode == 2, g += 5; end
  if inode == 3, g += 6; end
  g += 1;  % convert to 1-based
end

for iel = 0:nel-1
  for inode = 0:3
    for jnode = 0:3
      if locmatrix(inode+1, jnode+1) != 0
        A(conn(inode,iel), conn(jnode,iel)) += locmatrix(inode+1, jnode+1);
      end
    end
  end
end

% Penalty Dirichlet BCs at nodes 0 and 54 (1-based: 1 and 55)
A(1,  1)  += 1000;  b(1)  += 1000;
A(55, 55) += 1000;  b(55) += 1000;

x = A \ b;

printf ("x =\n");
printf ("%.15g\n", x);