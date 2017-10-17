function [ rate ] = r_dll_N2(Tin, xNH3,xNO, xO2 ) 
%UNTITLED Summary of this function goes here
%    T        NH3             NO   O2

[nn mm]=size(Tin);




% xNH3=max(1e-8,xNH3);
% xNO=max(1e-8,xNO);

H2Oin=0;

nline=length(Tin);
Xg(nn,6)=0;

%  NH3, O2, H2O, NO, N2, N2O, N2bal
for i=1:nline
Tin(i,1)=Tin(i);
Xg(i,1)=xNH3(i);
Xg(i,2)=xO2(i);
Xg(i,3)=H2Oin;
Xg(i,4)=xNO(i);
Xg(i,5)=0;
Xg(i,6)=0;
% Xg(i,7)=max(0.0, 1-sum(Xg(i,1:6)));
% Xg(i,7)=-1;
end

sfactors(1:18)=1;
%[rates,xsurf,elemrates] = nh3_ComputeRates(Tin,Xg,sfactors);
[rates,xsurf,elemrates] = chemistry(Tin,Xg);

%  NH3, O2, H2O, NO, N2, N2O, N2bal
% yNH3=rates(:,1);
yN2=rates(:,5);
% yN2O=rates(:,6);
% yNO=rates(:,4);
rate=yN2;

% fid2 = fopen('monitor_Chemistry.txt','a');
% for i=1:nn*mm
%  fprintf(fid2,'T= %g xNH3= %g xNO= %g xO2= %g rate=%g \n',Tin(i), xNH3(i),xNO(i), xO2(i),rate(i));
% end
% fclose(fid2);
   
    
rate=reshape(rate(:), nn,mm);

end

