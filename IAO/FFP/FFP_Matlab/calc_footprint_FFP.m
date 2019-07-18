function [FFP,flag_err]=calc_footprint_FFP(zm,z0,umean,h,ol,sigmav,ustar,varargin)

% [FFP,flag_err]=calc_footprint_FFP(zm,z0,umean,h,ol,sigmav,ustar,varargin)
% Derive a flux footprint estimate based on the simple parameterisation FFP
% 
% See Kljun, N., P. Calanca, M.W. Rotach, H.P. Schmid, 2015: 
% The simple two-dimensional parameterisation for Flux Footprint Predictions FFP.
% Geosci. Model Dev. 8, 3695-3713, doi:10.5194/gmd-8-3695-2015, for details.
% contact: n.kljun@swansea.ac.uk
%
%
% FFP Input
%    All inputs as scalars
%    zm        = Measurement height above displacement height (i.e. z-d) [m]
%    z0        = Roughness length [m] - enter [NaN] if not known 
%    umean     = Mean wind speed at zm [ms-1] - enter [NaN] if not known 
%                Either z0 or umean is required. If both are given,
%                z0 is selected to calculate the footprint
%    h         = Boundary layer height [m]
%    ol        = Obukhov length [m]
%    sigmav    = Standard deviation of lateral velocity fluctuations [ms-1]
%    ustar     = Friction velocity [ms-1]
%
%    Optional input (varargin):
%    Enter as calc_footprint_FFP(...,'OptionalInput',InputValue)
%    wind_dir  = Wind direction in degrees (of 360) for rotation of the footprint     
%    r         = Percentage of source area for which to provide contours, must be between 10% and 90%.
%                Can be either a single value (e.g., "80") or an array of percentage values 
%                (e.g., "[10:10:80]") 
%                Expressed either in percentages ("80") or in fractions of 1 ("0.8")
%                Default is [10:10:80]. Set to "NaN" for no output of percentages
%    nx        = Integer scalar defining the number of grid elements of the scaled footprint.
%                Large nx results in higher spatial resolution and higher computing time.
%                Default is 1000, nx must be >=600.
%    rslayer   = Calculate footprint even if zm within roughness sublayer: set rslayer = 1
%                Note that this only gives a rough estimate of the footprint as the model is not valid within 
%                the roughness sublayer. Default is 0 (i.e. no footprint for within RS).
%                z0 is needed for estimation of the RS.
%    crop      = Crop output area to size of the 80% footprint or the largest r given if crop=1
%
% FFP output
%    FFP          = Structure array with footprint data for measurement at [0 0 zm] m
%    FFP.x_ci_max = x location of footprint peak (distance from measurement) [m]
%    FFP.x_ci     = x array of crosswind integrated footprint [m]
%    FFP.f_ci     = array with footprint function values of crosswind integrated footprint [m-1] 
%    FFP.x_2d     = x-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
%    FFP.y_2d     = y-grid of 2-dimensional footprint [m], rotated if wind_dir is provided
%    FFP.f_2d     = Footprint function values of 2-dimensional footprint [m-2]
%    FFP.r        = Percentage of footprint as in input, if provided
%    FFP.fr       = Footprint value at r, if r is provided
%    FFP.xr       = x-array for contour line of r, if r is provided
%    FFP.yr       = y-array for contour line of r, if r is provided
%                   For array of percentage values, structure entries can be accessed 
%                   as FFP(1).r, FFP(1).xr, etc.
%    flag_err     = 1 in case of error, 0 otherwise
%
%
% Example
%    [FFP,flag_err]=calc_footprint_FFP(20,0.01,NaN,2000,-10,0.9,0.5,'wind_dir',30,'r',[10:20:80])
%
%
% created: 15 April 2015 natascha kljun
% version: 1.22
% last change: 29/09/2016 natascha kljun
%
% Copyright (C) 2015, Natascha Kljun
  
%--------------------------------------------------------------------
% Check input variables
%--------------------------------------------------------------------
opt_wind_dir = [];
opt_r        = [];
opt_nx       = [];
opt_rslayer  = [];
opt_crop     = [];

if nargin>=7 && nargin <=17
    if nargin > 7
        for f = 1:2:length(varargin)
            switch varargin{f}
                case 'wind_dir'
                    opt_wind_dir = varargin{f+1};
                case 'r'
                    opt_r = varargin{f+1};
                case 'nx'
                    opt_nx = varargin{f+1};
                case 'rslayer'
                    opt_rslayer = varargin{f+1};
                case 'crop'
                    opt_crop = varargin{f+1};
                otherwise
                    display(['Ignored unrecognised optional input "' varargin{f} '"'])
            end
        end
    end
    [ind_return,flag_err,wind_dir,r,nx,crop] = checkinput(zm,z0,umean,h,ol,sigmav,ustar, ...
         opt_wind_dir,opt_r,opt_nx,opt_rslayer,opt_crop);
else
    display('wrong number of input arguments')
    ind_return = 1;
    flag_err   = 1;
end

%--------------------------------------------------------------------
% Create output structure array
%--------------------------------------------------------------------
FFP = struct('x_ci_max',[],'x_ci',[],'f_ci',[],'x_2d',[],'y_2d',[],'f_2d',[], ...
             'r',[],'fr',[],'xr',[],'yr',[]);
if ind_return
    return
end %if

%--------------------------------------------------------------------
% Initialize model variables
%--------------------------------------------------------------------
a = 1.4524;
b = -1.9914;
c = 1.4622;
d = 0.1359;

ac = 2.17; 
bc = 1.66;
cc = 20.0;

xstar_end = 30;

%limit for neutral scaling
ol_n = 5000;

%von Karman
k = 0.4;

%--------------------------------------------------------------------
% Create scaled X* for crosswind integrated footprint
%--------------------------------------------------------------------
xstar_ci_param = linspace(d,xstar_end,nx+2);
xstar_ci_param = xstar_ci_param(2:end);
         
%--------------------------------------------------------------------
% Calculate crosswind integrated scaled F* 
%--------------------------------------------------------------------
fstar_ci_param = a.*(xstar_ci_param-d).^b .* exp(-c./(xstar_ci_param-d));     
ind_notnan     = ~isnan(fstar_ci_param);
fstar_ci_param = fstar_ci_param(ind_notnan);
xstar_ci_param = xstar_ci_param(ind_notnan);

%--------------------------------------------------------------------
% Calculate scaled sig_y*
%--------------------------------------------------------------------
sigystar_param = ac.*sqrt(bc.*(xstar_ci_param).^2 ./ (1+cc.*(xstar_ci_param)));

%--------------------------------------------------------------------
% Calculate real scale x and f_ci
%--------------------------------------------------------------------
if ~isnan(z0)
    if ol <=0 || ol >=ol_n
        xx  = (1 - 19.0.*zm./ol).^0.25;
        psi_f = log((1+xx.^2)./2) + 2.*log((1+xx)./2) - 2.*atan(xx) + pi./2;
    elseif ol > 0 && ol < ol_n
        psi_f = -5.3.*zm./ol;
    end
    
    x = xstar_ci_param.*zm ./ (1-(zm./h)) .* (log(zm./z0)-psi_f);
    if (log(zm./z0)-psi_f)>0
        x_ci = x;
        f_ci = fstar_ci_param./zm .* (1-(zm./h)) ./ (log(zm./z0)-psi_f);
    else
        flag_err = 1;
    end
else
    x = xstar_ci_param.*zm ./ (1-(zm./h)) .* (umean./ustar.*k);

    if (umean/ustar)>0
        x_ci = x;
        f_ci = fstar_ci_param./zm .* (1-(zm./h)) ./ (umean./ustar.*k);
    else
        flag_err = 1;
    end
end %if

if flag_err == 0
%--------------------------------------------------------------------
% Calculate maximum location of influence (peak location)
%--------------------------------------------------------------------
   xstarmax = -c./b+d;
   if ~isnan(umean)
       x_ci_max = xstarmax.*zm ./ (1-(zm./h)) .* (umean./ustar.*k);
   else
       x_ci_max = xstarmax.*zm ./ (1-(zm./h)) .* (log(zm./z0)-psi_f);
   end %if

%--------------------------------------------------------------------
% Calculate real scale sigy
%--------------------------------------------------------------------
   if abs(ol) >ol_n
       ol = -1E6;
   end
   if ol <=0 %convective
       scale_const = 1E-5.*abs(zm./ol).^(-1)+0.8;
   elseif ol > 0  %stable
       scale_const = 1E-5.*abs(zm./ol).^(-1)+0.55;
   end %if
   if scale_const>1
       scale_const = 1.0;
   end %if
   sigy         = sigystar_param./scale_const .*zm .*sigmav./ustar;
   sigy(sigy<0) = NaN;

%--------------------------------------------------------------------
% Calculate real scale f(x,y)
%--------------------------------------------------------------------
   dx    = x_ci(3)-x_ci(2);
   y_pos = 0:dx:(length(x_ci)/2)*dx*1.5;
   f_pos = NaN.*ones(length(f_ci),length(y_pos));
   for i=1:length(f_ci)
       f_pos(i,:) = f_ci(i) ./(sqrt(2.*pi).*sigy(i)) .* exp(-y_pos.^2./(2.*sigy(i).^2));
   end
   
%--------------------------------------------------------------------
% Complete footprint for negative y (symmetrical)
%--------------------------------------------------------------------
   y_neg = -fliplr(y_pos); 
   f_neg = fliplr(f_pos);
   y     = [y_neg(1:end-1) y_pos];
   f     = [f_neg(:,1:end-1) f_pos];

%--------------------------------------------------------------------
% Matrices for output
%--------------------------------------------------------------------
   [x_2d,y_2d] = meshgrid(x,y);
   f_2d = f';

%--------------------------------------------------------------------
% Derive footprint ellipsoid incorporating R% of the flux
% starting at peak value, if requested
%--------------------------------------------------------------------
   if ~isnan(r)
       rs = r;
   else
       if crop==1
           rs = 0.8;
       else
           rs = NaN;
       end
   end
   dy = dx;

   if ~isnan(rs)
       % Calculate integral of f_2d starting at peak value until R% are reached
       f_array = reshape(f_2d,1,size(f_2d,1)*size(f_2d,2));
       f_sort  = sort(f_array,'descend');
       f_sort  = f_sort(~isnan(f_sort));
       f_cum   = cumsum(f_sort).*dx.*dy;
       for i = 1:length(rs)
           f_diff     = abs(f_cum-rs(i));
           [~, ind_r] = min(f_diff);
           fr         = f_sort(ind_r);
           contour_r  = contourc(x,y,f_2d,[fr fr]);
    
           % Contourc adds info on level and number of vertices - replace with NaN
           ind_nan = find(contour_r(1,:)==fr);
           contour_r(:,ind_nan) = NaN;
           
           % Decrease number of digits and sort/unique
           try
               % matlab release 2014b or newer only
               contour_r = round(contour_r,1);
           catch
               contour_r = round(10.*contour_r)./10;
           end
       
           if ~isnan(r)               
               % Fill output structure
               FFP(i).r   = rs(i);
               FFP(i).fr  = fr;
               FFP(i).xr  = contour_r(1,:);
               FFP(i).yr  = contour_r(2,:);
           end
       end %for i
   end

%--------------------------------------------------------------------
% Crop domain
%--------------------------------------------------------------------
   if crop==1
       dminx = floor(min(contour_r(1,:)));
       dmaxx = ceil(max(contour_r(1,:)));
       dminy = floor(min(contour_r(2,:)));
       dmaxy = ceil(max(contour_r(2,:)));
       indy  = find(y_2d(:,1)>=dminy & y_2d(:,1)<=dmaxy);
       % extend by one row/column
       indy = [min(indy)-1; indy; max(indy)+1];
       indy = indy(indy>0 & indy<=size(x_2d,1));
       x_2d = x_2d(indy,:);
       y_2d = y_2d(indy,:);
       f_2d = f_2d(indy,:);
       indx = find(x_2d(1,:)>=dminx & x_2d(1,:)<=dmaxx);
       % extend by one row/column
       indx = [min(indx)-1 indx max(indx)+1];
       indx = indx(indx>0 & indx<=size(x_2d,2));
       x_2d = x_2d(:,indx);
       y_2d = y_2d(:,indx);
       f_2d = f_2d(:,indx);
   end

%--------------------------------------------------------------------
% Rotate footprint if requested
%--------------------------------------------------------------------
   if ~isnan(wind_dir)
       wind_dir_rad = wind_dir.*pi./180;
       dist         = sqrt(x_2d.^2 + y_2d.^2);
       angle        = atan2(y_2d, x_2d);
       x_2d_rot     = dist.*sin(wind_dir_rad - angle);
       y_2d_rot     = dist.*cos(wind_dir_rad - angle);
       
       if ~isnan(r)
           for i = 1:length(r)
               dist      = sqrt(FFP(i).xr.^2 + FFP(i).yr.^2);
               angle     = atan2(FFP(i).yr, FFP(i).xr);
               % Fill output structure
               FFP(i).xr = dist.*sin(wind_dir_rad - angle);
               FFP(i).yr = dist.*cos(wind_dir_rad - angle);
           end %for i
       end
   end %if

%--------------------------------------------------------------------
% Fill output structure
%--------------------------------------------------------------------
   FFP(1).x_ci_max = x_ci_max;
   FFP(1).x_ci     = x_ci;
   FFP(1).f_ci     = f_ci;
   if isnan(wind_dir)
       FFP(1).x_2d = x_2d;
       FFP(1).y_2d = y_2d;
       FFP(1).f_2d = f_2d;
   else
       FFP(1).x_2d = x_2d_rot;
       FFP(1).y_2d = y_2d_rot;
       FFP(1).f_2d = f_2d;
   end

end %if

end

%--------------------------------------------------------------------
% Function checkinput
%--------------------------------------------------------------------
function [ind_return,flag_err,wind_dir,r,nx,crop] = checkinput(zm,z0,umean, ...
          h,ol,sigmav,ustar,wind_dir,r,nx,rslayer,crop)

flag_err   = 0;
ind_return = 0;

if ~isscalar(ustar)
    display('Input of scalars only')
end

if zm<=0
    display('zm needs to be larger than 0 m')
    ind_return = 1;
elseif h<10
    display('h needs to be larger than 10 m')
    ind_return = 1;
elseif sigmav<0
    display('sig_v needs to be larger than 0 m/s')
    ind_return = 1;
elseif ustar<0
    display('ustar needs to be larger than m/s')
    ind_return = 1;
elseif zm>h
    display('zm needs to be smaller than h')
    ind_return = 1;
elseif zm/ol<=-15.5
    display('zm/L needs to be equal or larger than -15.5')
    ind_return = 1;
end
if isempty(wind_dir)
    wind_dir = NaN;
elseif wind_dir>360
    display('problem with wind direction')
    ind_return = 1;
elseif wind_dir<0
    display('problem with wind direction')
    ind_return = 1;
end

if isempty(r)
    r = 10:10:80;
end
if ~isnan(r(1))
    if max(r(:))>1
        r = r./100;
    end
    if max(r(:))>0.9
        display('R must be ,<= 0.9 or <=90%, larger values were removed')
        r = r(r<=0.9);
    end
    r = sort(r);
end

if isempty(nx)
    nx = 1000;
elseif nx<600
    display('nx must be >= 600')
    ind_return = 1;
end 
if isempty(rslayer)
    rslayer = 0;
end 

if ~isnan(z0)
    if z0 < 0
        display('z0 needs to be larger than 0')
        ind_return = 1;
    elseif zm<12.5*z0 && rslayer ~=1
        %changed to lowest limit of roughness sublayer definition
        display('zm needs to be above roughness sublayer')
        ind_return = 1;
    end
elseif ~isnan(umean)
    if umean < 0
        display('umean needs to be larger than 0')
        ind_return = 1;
    end
elseif isnan(z0) && isnan(umean)
    display('enter either z0 or umean')
    ind_return = 1;
end 

if ind_return
    flag_err = 1;
end

end

