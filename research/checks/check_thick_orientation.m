% Check: psm.m nZ line with a COLUMN thickness vector and >= 3 layers.
thick = [0.0733; 0.031; 0.054]; zOffset = 0.01; dzl = [1;2;3]*1e-4; nL = 3;
try
  nZ = ceil(sum([thick(1)-zOffset; thick(2:end).']./dzl(:)))+nL;
  printf('column thick: OK nZ=%d\n', nZ);
catch err
  printf('column thick: ERROR: %s\n', err.message);
end
thick = thick.';
nZ = ceil(sum([thick(1)-zOffset; thick(2:end).']./dzl(:)))+nL; printf('row thick: OK nZ=%d\n', nZ);
