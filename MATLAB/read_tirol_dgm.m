function [xx,yy,zz] = read_tirol_dgm(xlim,ylim,varargin)
% ==============================================================================
% Manuela Lehner
% March 2018
%
% INPUT: map limits in MGI coordinates (Austria GK West)
%        xlim, ylim ... 2 component vectors
%        OPTIONAL: 'dom' ... read DOM instead of DGM
%        OPTIONAL: '10m' ... read 10 m data istead of 5 m
% OUTPUT: xx,yy ... coordinate vectors
%         zz   ... elevation array
% ==============================================================================
n = 4;
f = waitbar(1/n,"Loading data");
% -----
if nargin<2; error('Not enough input arguments'); end
infile = 'DGM_Tirol_5m_epsg31254';
for iv = 1 : length(varargin)
  if strcmp(varargin{iv}, 'dom')
    infile = 'DOM_Tirol_5m_epsg31254';
  elseif strcmp(varargin{iv}, '10m')
    infile = '~/topo/Tirol_DOM5/DGM_Tirol_10m_epsg31254';
  end
end

% ----- read file
f = waitbar(2/n,f,"Reading file");

[zz, R] = geotiffread(infile);

f = waitbar(3/n,f,"Calculating coordinates");

% ----- calculate coordinates
xx = R.XWorldLimits(1) + [1:R.RasterSize(2)].*R.CellExtentInWorldX;
yy = fliplr(R.YWorldLimits(1) + [1:R.RasterSize(1)].*R.CellExtentInWorldY);
if( abs(xx(end)-R.XWorldLimits(2))>1.0 | abs(yy(1)-R.YWorldLimits(2))>1.0 )
  error('Coordinate calculation');
end
f = waitbar(4/n,f,"Extracting domain");

% ----- extract domain specified by xlim and ylim
indx = find(xx>=xlim(1) & xx<=xlim(2));
indy = find(yy>=ylim(1) & yy<=ylim(2));
xx = xx(indx);
yy = yy(indy);
zz = zz(indy,indx);
close(f)