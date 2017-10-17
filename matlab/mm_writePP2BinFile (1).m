function [ ] = mm_writePP2BinFile( zpp, name )

coefs(:,1) = zpp.coefs(:);

% Coefs
h = [name, '_coefs.bin'];
fid = fopen(h,'wb');
fwrite(fid,length(coefs),'real*8'); % This works with ifort, although an integer is expected. Other compilers will have problems with that.
fwrite(fid,coefs,'real*8');
fclose(fid);
