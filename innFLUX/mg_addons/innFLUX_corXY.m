TRACER_x = 17;  % index of x-VOC (note in ptrmsYYYYMMDD.mat EC.data(:,1) is time)
TRACER_y = 27;  % index of y-VOC (note in ptrmsYYYYMMDD.mat EC.data(:,1) is time)

folder_paths = {'G:\IAO2018_PTR\IAO-preIOP2018\NO_mode\vC\concentrations';...
    'G:\IAO2018_PTR\IAO-IOP2018_Summer\NO_mode\vC\concentrations'};
prefixes{1} = 'ptrms';
suffix = '.mat';
date_format = 'yyyymmdd';

% plotting individual 1/2hrs for debugging
DoPlot = false;
if DoPlot
    figure('Position',[50 50 1800 600])
end
% smoothing before regression
boxwidth = 101; % must be odd integer

files = list_files_datestring_paths(folder_paths, prefixes, suffix, date_format);
EC = load([files(1).path '\' files(1).name]);

nums_T = size(EC.data,2)-1;
names_T = EC.header;

corXY = NaN(nums_T,nums_T,length(files)*48); % 3D: will hold the correlation matrix (corrcoef)
p_XY = corXY;
corXY_LO = corXY;
corXY_UP = corXY;
utc = NaN(length(files)*48,1);

%IsFirst = true;

for f_dex = 1:length(files)
    EC = load([files(f_dex).path '\' files(f_dex).name]);
   
        % get complete time series of m92 to check spikes at file edges
        if(0)
        if f_dex == 1
            TS.time = EC.data(:,1);
            TS.m93 = EC.data(:,53);
        else
            TS.time = cat(1,TS.time,EC.data(:,1));
            TS.m93 = cat(1,TS.m93,EC.data(:,53));
        end
if(1)
    for i=1:48;
        idx=(1:18000)+(i-1)*18000;
        % selecting later % idx(~isfinite(EC.data(idx,TRACER_x))|~isfinite(EC.data(idx,TRACER_y)))=[];
        
        
        
        % smooth and trim
        Tmat = EC.data(idx,2:end);
        Ttime = EC.data(idx,1);
        

        
        if(boxwidth>1)
            %width = oddify(boxwidth);
            %   fast_aux = smooth(fast_sig,width);  replace smoothing by filter
            Tmat = filter(ones(boxwidth,1)/boxwidth,1,Tmat);
            Tx(1:boxwidth,:) = NaN;    %   filter needs 'width' datapoints to spin up
%             Ty = filter(ones(boxwidth,1)/boxwidth,1,Ty);
%             Ty(1:boxwidth) = NaN;    %   filter needs 'width' datapoints to spin up
        end
        
        % trim first four and last minute of each 1/2hr
        Tmat([1:4*600 18000-600:18000],:) = NaN;
        %Ty([1:600 18000-600:18000]) = NaN;
        
        T_nan = find(any(~isfinite(Tmat),2));
        Tmat(T_nan,:)=[];   % remove rows where there is any NaN
        
        if numel(Tmat)>1
            [cc,cp,lo,up] = corrcoef(Tmat);
            corXY(1:nums_T,1:nums_T,(i+(f_dex-1)*48)) = cc;
            p_XY(1:nums_T,1:nums_T,(i+(f_dex-1)*48)) = cp;
            corXY_LO(1:nums_T,1:nums_T,(i+(f_dex-1)*48)) = lo;
            corXY_UP(1:nums_T,1:nums_T,(i+(f_dex-1)*48)) = up;

            
            utc(i+(f_dex-1)*48)=min(Ttime);
            
            if DoPlot
                
                subplot(1,3,1:2)
                plot(Ttime(T_dex),Tx(T_dex),'r')
                hold on
                plot(Ttime(T_dex),Ty(T_dex),'k')
                hold off
                xlim([Ttime(1),Ttime(end)])
                datetick('x','keeplimits')
                grid on
                xlabel('UTC','FontSize',14)
                ylabel('VMR [ppb]','FontSize',14)
                legend(S_leg)
                
                subplot(1,3,3)
                plot(Tx(T_dex),Ty(T_dex),'.')
                grid on
                xlabel(name_x,'FontSize',14)
                ylabel(name_y,'FontSize',14)
                S_reg = sprintf('y = %0.2f * x + %0.2f; r = %0.2f',reg.b,reg.a,reg.r);
                legend(S_reg,'FontSize',14,'Location','NW')
            end
            
        end
    end
end
end

if(1)
r = squeeze(corXY(TRACER_x-1,TRACER_y-1,:));
name_x = char(names_T(TRACER_x));
name_y = char(names_T(TRACER_y));

figure('Position', [50 50 2400 900])
subplot(1,2,1)

plot(utc,r,'.')
datetick('x','keeplimits')
xlabel('UTC','FontSize',14)
ylabel('correlation','FontSize',14)
title(['Hourly Regression Coefficient(' char(name_x) ',' char(name_y) ')'],'FontSize',16)
subplot(1,2,2)
boxplot(r,ceil(mod(utc,1)*24))

end