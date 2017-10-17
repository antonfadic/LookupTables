function [ ] = mm_writeRISTO2BinFile( LookUpz, varargin )
% This writes the fungrid_* and dergrid_* files

Options.writeFunGrid = 1;
Options.writeDerGrid = 1;

for i = 1:2:length(varargin)
    Options.(varargin{i}) = varargin{i+1};
end

for luIndex = 1:length(LookUpz);
    % fungrid schreiben
    if ( Options.writeFunGrid == 1 )
        h = ['fungrid_',LookUpz(luIndex).Name,'.bin'];
        fid = fopen(h,'wb');
        %fwrite(fid, length(LookUpz(luIndex).OferSpline.fungrid), 'real*8');
        fwrite(fid, LookUpz(luIndex).OferSpline.fungrid, 'real*8');
        fclose(fid);
    end
    
    % dergrid schreiben
    if ( Options.writeDerGrid == 1 )
        h = ['dergrid_',LookUpz(luIndex).Name,'.bin'];
        fid = fopen(h,'w');
        %fwrite(fid ,length(LookUpz(luIndex).OferSpline.dergrid), 'real*8');
        fwrite(fid, LookUpz(luIndex).OferSpline.dergrid, 'real*8');
        fclose(fid);
    end
end
