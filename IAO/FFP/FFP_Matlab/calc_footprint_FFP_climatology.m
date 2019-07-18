function [FFP,flag_err]=calc_footprint_FFP_climatology(zm,z0,umean,h,ol,sigmav,ustar,wind_dir,varargin)

% [FFP,flag_err]=calc_footprint_FFP_climatology(zm,z0,umean,h,ol,sigmav,ustar,wind_dir,varargin)
% Derive a flux footprint climatology based on the simple parameterisation FFP
% 
% See Kljun, N., P. Calanca, M.W. Rotach, H.P. Schmid, 2015: 
% The simple two-dimensional parameterisation for Flux Footprint Predictions FFP.
% Geosci. Model Dev. 8, 3695-3713, doi:10.5194/gmd-8-3695-2015, for details.
% contact: n.kljun@swansea.ac.uk
%
% This function calculates footprints within a fixed physical domain for a series of
% time steps, rotates footprints into the corresponding wind direction and aggregates 
% all footprints to a footprint climatology. The percentage of source area is
% calculated for the footprint climatology.
% For determining the optimal extent of the domain (large enough to include footprints)
% use calc_footprint_FFP.m
%
% FFP Input
%    All vectors need to be of equal length (one value for each time step)
%    zm       = Measurement height above displacement height (i.e. z-d) [m]
%               usually a scalar, but can also be a vector 
%    z0       = Roughness length [m] - enter [NaN] if not known 
%               usually a scalar, but can also be a vector 
%    umean    = Vector of mean wind speed at zm [ms-1] - enter [NaN] if not known 
%               Either z0 or umean is required. If both are given, z0 is selected to calculate the footprint
%    h        = Vector of boundary layer height [m]
%    ol       = Vector of Obukhov length [m]
%    sigmav   = Vector of standard deviation of lateral velocity fluctuations [ms-1]
%    ustar    = Vector of friction velocity [ms-1]
%    wind_dir = Vector of wind direction in degrees (of 360) for rotation of the footprint     
%
%    Optional input (varargin):
%    Enter as calc_footprint_FFP_climatology(...,'OptionalInput',InputValue)
%    domain       = Domain size as an array of [xmin xmax ymin ymax] [m]
%                   Footprint will be calculated for a measurement at [0 0 zm] m
%                   Default is smallest area including the r% footprint or [-1000 1000 -1000 1000]m,
%                   which ever smallest (80% footprint if r not given)
%    dx, dy       = Cell size of domain [m]
%                   Small dx,dy result in higher spatial resolution and higher computing time
%                   Default is dx = dy = 2 m (if neither domain nor nx and ny are given). 
%                   If only dx is given, dx=dy.
%    nx, ny       = Two integer scalars defining the number of grid elements in x and y
%                   Large nx and ny result in higher spatial resolution and higher computing time
%                   Default is nx = ny = 1000. If only nx is given, nx=ny 
%                   If dx,dy and nx,ny are given, dx,dy is given priority
%    r            = Percentage of source area for which to provide contours, must be between 10% and 90%.
%                   Can be either a single value (e.g., "80") or an array of percentage values (e.g., "[10:10:80]") 
%                   Expressed either in percentages ("80") or in fractions of 1 ("0.8")
%                   Default is [10:10:80]. Set to "NaN" for no output of percentages
%    rslayer      = Calculate footprint even if zm within roughness sublayer: set rslayer = 1
%                   Note that this only gives a rough estimate of the footprint as the model is not valid within 
%                   the roughness sublayer. Default is 0 (i.e. no footprint for within RS).
%                   z0 is needed for estimation of the RS.
%    smooth_data  = Apply convolution filter to smooth footprint climatology if smooth_data=1 (default)
%    crop         = Crop output area to size of the 80% footprint or the largest r given if crop=1
%    pulse        = Display progress of footprint calculations every pulse-th footprint (e.g., "100")
%
% FFP output
%    FFP          = Structure array with footprint climatology data for measurement at [0 0 zm] m
%    FFP.x_2d     = x-grid of footprint climatology [m]
%    FFP.y_2d     = y-grid of footprint climatology [m]
%    FFP.fclim_2d = Normalised footprint function values of footprint climatology [m-2]
%    FFP.r        = Percentage of footprint as in input, if provided
%    FFP.fr       = Footprint value at r, if r is provided [m-2]
%    FFP.xr       = x-array for contour line of r, if r is provided [m]
%    FFP.yr       = y-array for contour line of r, if r is provided [m]
%                   For array of percentage values, structure entries can be accessed 
%                   as FFP(1).r, FFP(1).xr, etc.
%    FFP.n        = Number of footprints calculated and included in footprint climatology
%    flag_err     = 0 if no error, 1 in case of error, 2 if not all contour plots (r%) within specified domain
%                 
%                   
% Example
%    zm=20; z0=0.01; umean=NaN;
%    h=[2000 1800 1500]; ol=[-10 -100 -500]; sigmav=[0.9 0.7 0.3]; ustar=[0.5 0.3 0.4]; wind_dir=[30 50 70];
%    [FFP,flag_err] = calc_footprint_FFP_climatology(zm,z0,umean,h,ol,sigmav,ustar,wind_dir, ...
%                     'domain',[-100 1000 -100 1000],'nx',1100,'r',[10:10:80],'smooth_data',1)
%    
%
% created: 19 February 2016 natascha kljun
% version: 1.22
% last change: 18/09/2016 natascha kljun
%
% Copyright (C) 2015, Natascha Kljun

%--------------------------------------------------------------------
% Check input variables
%--------------------------------------------------------------------
opt_domain  = [];
opt_dx      = [];
opt_dy      = [];
opt_nx      = [];
opt_ny      = [];
opt_r       = [];
opt_rslayer = [];
opt_smooth  = [];
opt_crop    = [];
opt_pulse   = [];

if nargin >= 8 && nargin<=28
    if nargin > 8
        for f = 1:2:length(varargin)
            switch varargin{f}
                case 'domain'
                    opt_domain = varargin{f+1};
                case 'dx'
                    opt_dx = varargin{f+1};
                case 'dy'
                    opt_dy = varargin{f+1};
                case 'nx'
                    opt_nx = varargin{f+1};
                case 'ny'
                    opt_ny = varargin{f+1};
                case 'r'
                    opt_r = varargin{f+1};
                case 'rslayer'
                    opt_rslayer = varargin{f+1};
                 case 'smooth_data'
                    opt_smooth = varargin{f+1};
                case 'crop'
                    opt_crop = varargin{f+1};
                case 'pulse'
                    opt_pulse = varargin{f+1};
                otherwise
                    display(['Ignored unrecognised optional input "' varargin{f} '"'])
            end
        end
    end
    [ind_return,flag_err,valid,ts_len,zm,z0,wind_dir,xmin,xmax,ymin,ymax, ...
        nx,ny,dx,dy,r,smooth_data,crop,pulse] = checkinput(zm,z0,umean,h,ol,sigmav,ustar, ...
        wind_dir,opt_domain,opt_dx,opt_dy,opt_nx,opt_ny,opt_r,opt_rslayer,opt_smooth,opt_crop,opt_pulse);
else
    display('Wrong number of input arguments')
    ind_return = 1;
    flag_err   = 1;
end

%--------------------------------------------------------------------
% Create output structure array
%--------------------------------------------------------------------
FFP = struct('x_2d',[],'y_2d',[],'fclim_2d',[],'r',[],'fr',[],'xr',[],'yr',[],'n',[]);
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

%limit for neutral scaling
ol_n = 5000;

%von Karman
k = 0.4;

%--------------------------------------------------------------------
% Define domain
%--------------------------------------------------------------------
% Define physical domain in cartesian and polar coordinates
% Cartesian coordinates
x           = linspace(xmin,xmax,nx+1);
y           = linspace(ymin,ymax,ny+1);
[x_2d,y_2d] = meshgrid(x,y);
% Polar coordinates
% Set theta such that North is pointing upwards and angles increase clockwise
rho   = sqrt(x_2d.^2 + y_2d.^2);
theta = atan2(x_2d,y_2d);

% initialize raster for footprint climatology
fclim_2d = zeros(size(x_2d));

%--------------------------------------------------------------------
% Start loop on time series
%--------------------------------------------------------------------
for foot_loop = 1:ts_len
    
    if mod(foot_loop,pulse) == 0
        display(['Calculating footprint ' int2str(foot_loop) ' of ' int2str(ts_len)])
    end %if
    
    if valid(foot_loop)
        
%--------------------------------------------------------------------
% Create local (within loop) variables        
%--------------------------------------------------------------------
        wind_dirl = wind_dir(foot_loop);
        oll       = ol(foot_loop);
        zml       = zm(foot_loop);
        if ~isnan(z0(foot_loop))
            z0l       = z0(foot_loop);
        elseif ~isnan(umean(foot_loop))
            umeanl    = umean(foot_loop);
        end
        hl        = h(foot_loop);
        ustarl    = ustar(foot_loop);
        sigmavl   = sigmav(foot_loop);
        
%--------------------------------------------------------------------
% Rotate coordinates into wind direction
%--------------------------------------------------------------------
        wind_dir_rad = wind_dirl.*pi./180;
        thetal       = theta-wind_dir_rad;
 
%--------------------------------------------------------------------
% Create real scale crosswind integrated footprint and dummy for
% rotated scaled footprint
%--------------------------------------------------------------------
        fstar_ci_dummy = zeros(size(x_2d));
        f_ci_dummy     = zeros(size(x_2d));
        if min(z0)>0
            if oll <= 0 || oll >= ol_n
                xx    = (1 - 19.0.*zml./oll).^0.25;
                psi_f = log((1+xx.^2)./2) + 2.*log((1+xx)./2) - 2.*atan(xx) + pi./2;
            elseif oll > 0 && oll < ol_n
                psi_f = -5.3.*zml./oll;
            end
            if (log(zml./z0l)-psi_f)>0
                xstar_ci_dummy     = rho.*cos(thetal)./zml .* (1-(zml./hl)) ./ (log(zml./z0l)-psi_f);
                px                 = xstar_ci_dummy>d;
                fstar_ci_dummy(px) = a.*(xstar_ci_dummy(px)-d).^b .* exp(-c./(xstar_ci_dummy(px)-d));
                f_ci_dummy(px)     = fstar_ci_dummy(px)./zml .* (1-(zml./hl)) ./ (log(zml./z0l)-psi_f);
            else
                flag_err = 1;
            end
        else
            xstar_ci_dummy     = rho.*cos(thetal)./zml .* (1-(zml./hl)) ./ (umeanl./ustarl.*k);
            px                 = xstar_ci_dummy>d;
            fstar_ci_dummy(px) = a.*(xstar_ci_dummy(px)-d).^b .* exp(-c./(xstar_ci_dummy(px)-d));
            f_ci_dummy(px)     = fstar_ci_dummy(px)./zml .* (1-(zml./hl)) ./ (umeanl./ustarl.*k);
        end %if
       
%--------------------------------------------------------------------
% Calculate dummy for scaled sig_y* and real scale sig_y
%--------------------------------------------------------------------
        sigystar_dummy     = zeros(size(x_2d));
        sigystar_dummy(px) = ac.*sqrt(bc.*(xstar_ci_dummy(px)).^2 ./ (1+cc.*(xstar_ci_dummy(px))));

        if abs(oll) >ol_n
            oll = -1E6;
        end
        if oll <=0 %convective
            scale_const = 1E-5.*abs(zml./oll).^(-1)+0.8;
        elseif oll > 0  %stable
            scale_const = 1E-5.*abs(zml./oll).^(-1)+0.55;
        end %if
        if scale_const>1
            scale_const = 1.0;
        end %if
        sigy_dummy     = zeros(size(x_2d));
        sigy_dummy(px) = sigystar_dummy(px)./scale_const .*zml .*sigmavl./ustarl;
        sigy_dummy(sigy_dummy<0) = NaN;

%--------------------------------------------------------------------
% Calculate real scale f(x,y)
%--------------------------------------------------------------------
        f_2d     = zeros(size(x_2d));
        f_2d(px) = f_ci_dummy(px) ./ (sqrt(2.*pi).*sigy_dummy(px)) .* ...
                   exp(-(rho(px).*sin(thetal(px))).^2./(2.*sigy_dummy(px).^2));
   
%--------------------------------------------------------------------
% Add to footprint climatology raster
%--------------------------------------------------------------------
        fclim_2d = fclim_2d + f_2d;

    end %if valid
end %for foot_loop

%--------------------------------------------------------------------
% Normalize and smooth footprint climatology
%--------------------------------------------------------------------
fclim_2d = fclim_2d./sum(valid);

if smooth_data
    skernel  = [0.05 0.1 0.05; 0.1 0.4 0.1; 0.05 0.1 0.05];
    fclim_2d = conv2(fclim_2d,skernel,'same');
    fclim_2d = conv2(fclim_2d,skernel,'same');
end %if

%--------------------------------------------------------------------
% Derive footprint ellipsoid incorporating R% of the flux, if requested,
% starting at peak value
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

if ~isnan(rs)
    % Calculate integral of fclim_2d starting at peak value until R% are reached
    f_array = reshape(fclim_2d,1,size(fclim_2d,1)*size(fclim_2d,2));
    f_sort  = sort(f_array,'descend');
    f_sort  = f_sort(~isnan(f_sort));
    f_cum   = cumsum(f_sort).*dx.*dy;
    for i = 1:length(rs)
        f_diff     = abs(f_cum-rs(i));
        [~, ind_r] = min(f_diff);
        fr         = f_sort(ind_r);
        contour_r  = contourc(x,y,fclim_2d,[fr fr]);
    
        % Contourc adds info on level and number of vertices - replace with NaN
        ind_nan = find(contour_r(1,:)==fr);
        contour_r(:,ind_nan) = NaN;

        % Decrease number of digits and sort/unique, matlab release 2014b or newer!
        try
            % matlab release 2014b or newer only
            contour_r = round(contour_r,1);
        catch
            contour_r = round(10.*contour_r)./10;
        end
    
        if ~isnan(r)
            % No output of R% contour lines if outside domain extent
            if max(contour_r(2,:))>=max(y_2d(:)) || max(contour_r(1,:))>=max(x_2d(:)) || ...
                    min(contour_r(2,:))<=min(y_2d(:)) || min(contour_r(1,:))<=min(x_2d(:))
                flag_err         = 2;
                contour_r        = [];
                contour_r(1:2,1) = NaN;
                fr               = NaN;
            end
    
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
    if isnan(contour_r(1,:))
        while isnan(contour_r(1,:))
            i = i-1;
            contour_r(1,1:length(FFP(i).xr)) = FFP(i).xr;
            contour_r(2,1:length(FFP(i).yr)) = FFP(i).yr;
        end
    end
    dminx = floor(min(contour_r(1,:)));
    dmaxx = ceil(max(contour_r(1,:)));
    dminy = floor(min(contour_r(2,:)));
    dmaxy = ceil(max(contour_r(2,:)));
    if dminy>=ymin && dmaxy<=ymax
        indy = find(y_2d(:,1)>=dminy & y_2d(:,1)<=dmaxy);
        % extend by one row/column
        indy = [min(indy)-1; indy; max(indy)+1];
        indy = indy(indy>0 & indy<=size(y_2d,1));
        x_2d     = x_2d(indy,:);
        y_2d     = y_2d(indy,:);
        fclim_2d = fclim_2d(indy,:);
    end
    if dminx>=xmin && dmaxx<=xmax
        indx = find(x_2d(1,:)>=dminx & x_2d(1,:)<=dmaxx);
        % extend by one row/column
        indx = [min(indx)-1 indx max(indx)+1];        
        indx = indx(indx>0 & indx<=size(x_2d,2));
        x_2d     = x_2d(:,indx);
        y_2d     = y_2d(:,indx);
        fclim_2d = fclim_2d(:,indx);
    end
end

%--------------------------------------------------------------------
% Fill output structure
%--------------------------------------------------------------------
FFP(1).x_2d     = x_2d;
FFP(1).y_2d     = y_2d;
FFP(1).fclim_2d = fclim_2d;
FFP(1).n        = sum(valid);

end
  

%--------------------------------------------------------------------
% Function checkinput
%--------------------------------------------------------------------
function [ind_return,flag_err,valid,ts_len,zm,z0,wind_dir,xmin,xmax,ymin,ymax, ...
         nx,ny,dx,dy,r,smooth_data,crop,pulse] = checkinput(zm,z0,umean,h,ol,sigmav,ustar, ...
         wind_dir,domain,dx,dy,nx,ny,r,rslayer,smooth_data,crop,pulse)
      
flag_err   = 0;
ind_return = 0;
ts_len     = length(ustar);

if numel(ustar) ~= ts_len
    display('Input of scalars or vectors only')
end
if isscalar(zm)
    zm = zm .*ones(size(ustar));
elseif length(zm) ~= ts_len 
   display('zm must be either scalar or of same length as other vectors')
   ind_return = 1;
end
if isscalar(z0)
    z0 = z0 .*ones(size(ustar));
elseif length(z0) ~= ts_len 
   display('z0 must be either scalar or of same length as other vectors')
   ind_return = 1;
end
if length(h) ~= ts_len || length(ol) ~= ts_len || length(sigmav) ~= ts_len || ... 
   length(wind_dir) ~= ts_len 
   display('Input vectors must be of same length')
   ind_return = 1;
end

if min(zm)<=0
    display('zm must be larger than 0')
    ind_return = 1;
elseif min(h)<10
    display('h must be larger than 10 m')
    ind_return = 1;
elseif min(sigmav)<0
    display('sig_v must be larger than 0')
    ind_return = 1;
elseif min(ustar)<0
    display('ustar must be larger than 0')
    ind_return = 1;
elseif max(wind_dir)>360
    display('problem with wind direction')
    ind_return = 1;
elseif min(wind_dir)<0
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

if isempty(rslayer)
    rslayer = 0;
end 

if isempty(smooth_data)
    smooth_data = 1;
elseif smooth_data ~= 0 && smooth_data ~= 1
    display('smooth_data must be 0 or 1')
    ind_return = 1;
end 

if isempty(crop)
    crop = 0;
end 

if isempty(pulse)
    if ts_len <=20
        pulse = 1;
    else
        pulse = round(ts_len/100);
    end
end 

if isempty(domain) || isnan(domain(1))
    xmin = -1000;
    xmax = 1000;
    ymin = -1000;
    ymax = 1000;
else
    if numel(domain)~=4
        display('domain must be an array of four elements [xmin xmax ymin ymax]')
        ind_return = 1;
    else
        min_extent = min(domain(4)-domain(3),domain(2)-domain(1));
        if min_extent<1
            display('domain extent must be larger than 1 m in both x and y')
            ind_return = 1;
        else
            xmin = round(domain(1));
            xmax = round(domain(2));
            ymin = round(domain(3));
            ymax = round(domain(4));
        end
    end
end %if
if ~isempty(dx) || ~isempty(dy)
    if isempty(dy)
        dy = dx;
    elseif isempty(dx)
        dx = dy;   
    end
    nx = round((xmax-xmin)./dx);
    ny = round((ymax-ymin)./dy);
    if max([nx ny])>2000
        display('very small dxy - may slow down calculation and cause problems when plotting')
    end
else
    if isempty(nx) && isempty(ny)
        nx = 1000;
        ny = nx;
    elseif isempty(ny)
        ny = nx;
    elseif isempty(nx)
        nx = ny;
    else
        if mod(nx,1) ~= 0 || mod(ny,1) ~= 0
            display('nx (and/or ny) are rounded to next integer values')
            nx = round(nx);
            ny = round(ny);
        elseif min(nx)<10 || min(ny)<10
            display('nx and ny extent must be larger than 10 m in both x and y')
            ind_return = 1;
            nx         = NaN;
            ny         = NaN;
        elseif max(nx)>2000 || max(ny)>2000 
            display('very large nx or ny - may slow down calculation and cause problems when plotting')
        end
    end
    dx = (xmax-xmin)./nx;
    dy = (ymax-ymin)./ny;
end

if sum(zm>h)>0
    display('zm must be smaller than h')
    ind_return = 1;
end

if ~isnan(z0)
    if min(z0) < 0
        display('z0 must be larger than 0')
        ind_return = 1;
    end
elseif ~isnan(umean)
    if min(umean) < 0
        display('umean must be larger than 0')
        ind_return = 1;
    end
elseif isnan(z0) && isnan(umean)
    display('enter either z0 or umean')
    ind_return = 1;
end

valid          = ones(size(ustar));
ind_nan        = isnan(ustar);
valid(ind_nan) = 0;
ind_nan        = ustar<0.1;
valid(ind_nan) = 0;
ind_nan        = isnan(ol);
valid(ind_nan) = 0;
ind_nan        = isnan(sigmav);
valid(ind_nan) = 0;
ind_nan        = isnan(wind_dir);
valid(ind_nan) = 0;

%remove single data points if following conditions not fulfilled

%zm within roughness sublayer (changed to lowest limit)
%can only check if z0 given - condition may be ignored if rs_layer=1
if ~isnan(z0) 
    ind_nan = (zm <= z0.*12.5);
    if rslayer~=1 && sum(ind_nan)>0
        valid(ind_nan)  = 0;
        display('zm must be above roughness sublayer')
    end
end
ind_nan = zm./ol<= -15.5;
if sum(ind_nan)>0
    valid(ind_nan) = 0;
    display('zm/L must be equal or larger than -15.5')
end

if sum(valid)<1
    ind_return = 1;
end

if ind_return
    flag_err = 1;
end 

end
