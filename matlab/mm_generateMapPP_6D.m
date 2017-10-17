function [ LookUpz ] = mm_generateMapPP_6D(polyorder, scaleType, varargin)

% The dimension are:
% Number    Name    lb      ub      scale 
% 1         T       300     800     reciprocal
% 2         CH4     1e-8    1e-3    logarithmic
% 3         H2O      1e-8    1e-3    logarithmic
% 4         CO      0.02    0.2     logarithmic
% 5         H2      0.02    0.2     logarithmic
% 6         CO2      0.02    0.2     logarithmic

if length(varargin) ~=0
    filename = varargin{1};
else
    filename =  'spline_4D_20_mal_20_Alex_minimal.mat';
end

% yy* contain the function values on the grid which is specified by r*.
load(filename);
% yyCH4 = -yyCH4;

nDims = 4;

%% Generate Spline maps - minimal version:
    %CH4 H2O CO H2 CO2
LookUpz(1).Name = 'NH3';
LookUpz(2).Name = 'N2';
LookUpz(3).Name = 'N2O';
% LookUpz(2).Name = 'H2O';
% LookUpz(3).Name = 'CO';
% LookUpz(4).Name = 'H2';
% LookUpz(2).Name = 'CO2';

k_thingy = cell(size(ones(1,nDims)));
for dimIndex = 1:nDims
    k_thingy{dimIndex} = polyorder;
end

for i = 1:length(LookUpz)
    yyvar = eval(['yy' LookUpz(i).Name]);
    switch scaleType
        case 'log'
            tempyy = log(abs(yyvar));
        case 'nolog'
            tempyy = yyvar;
    end

LookUpz(i).pp = fn2fm(spapi(k_thingy, {r1 r2 r3 r4}, tempyy), 'pp');
end
