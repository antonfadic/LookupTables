%% Load the grid file from C
clear all; close all; clc;

fileID = fopen('grid.bin');
A = fread(fileID,'double');
A = reshape(A,4,[]);
fclose(fileID);
%A(1,:) = 1./A(1,:);
%A(2:end,:) = exp(A(2:end,:));

%% Obtain the rates using A
%A(1,:) T
%A(2,:) yNH3
%A(3,:) yNO
%A(4,:) yO2
nDimOut=3;
LengthVal = length(A(1,:));
B = zeros(LengthVal,nDimOut);
clc;
tic
%call the mechanism
    for i=1:length(A(1,:))
        B(i,1) = r_dll_NH3(A(1,i),A(2,i),A(3,i),A(4,i));
        B(i,2) = r_dll_N2(A(1,i),A(2,i),A(3,i),A(4,i));
        B(i,3) = r_dll_N2O(A(1,i),A(2,i),A(3,i),A(4,i));
        if(mod(i,1000)==1)
            fprintf('%i \n',i);
        end
    end
    time=toc;

%% Save the grid file in binary format and load it later into C
fileID = fopen('Table.bin','wb');
B = reshape(B,1,[]);
surfDens= 2.7063e-6;
fwrite(fileID,log(B.*surfDens),'double');
fclose(fileID)