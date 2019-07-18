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

r = NaN(length(files)*48,1);
k = r;
d = r;
utc = r;

IsFirst = true;

for f_dex = 1:length(files)
    EC = load([files(f_dex).path '\' files(f_dex).name]);
    if IsFirst
        name_x = char(EC.header(TRACER_x));
        name_y = char(EC.header(TRACER_y));
        S_leg{1} = name_x;
        S_leg{2} = name_y;
        IsFirst = false;
    end
    
    for i=1:48;
        idx=(1:18000)+(i-1)*18000;
        % selecting later % idx(~isfinite(EC.data(idx,TRACER_x))|~isfinite(EC.data(idx,TRACER_y)))=[];
        
        
        
        % smooth and trim
        Tx = EC.data(idx,TRACER_x);
        Ty = EC.data(idx,TRACER_y);
        Ttime = EC.data(idx,1);
        
        if(boxwidth>1)
            %width = oddify(boxwidth);
            %   fast_aux = smooth(fast_sig,width);  replace smoothing by filter
            Tx = filter(ones(boxwidth,1)/boxwidth,1,Tx);
            Tx(1:boxwidth) = NaN;    %   filter needs 'width' datapoints to spin up
            Ty = filter(ones(boxwidth,1)/boxwidth,1,Ty);
            Ty(1:boxwidth) = NaN;    %   filter needs 'width' datapoints to spin up
        end
        
        % trim first and last minute of each 1/2hr
        Tx([1:600 18000-600:18000]) = NaN;
        Ty([1:600 18000-600:18000]) = NaN;
        T_dex = find(isfinite(Tx) & isfinite(Ty));
        
        if numel(T_dex)>0
            [c1,c2,c3] = regression(Tx(T_dex),Ty(T_dex),'one');
            k(i+(f_dex-1)*48)=c2;
            d(i+(f_dex-1)*48)=c3;
            r(i+(f_dex-1)*48)=c1;            
%             reg = reg_mg(Tx(T_dex),Ty(T_dex));
%             k(i+(f_dex-1)*24)=reg.b;
%             d(i+(f_dex-1)*24)=reg.a;
%             r(i+(f_dex-1)*24)=reg.r;
            utc(i+(f_dex-1)*48)=min(EC.data(idx,1));
            
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
figure('Position', [50 50 2400 900])
subplot(1,2,1)
plot(utc,r,'.')
datetick('x','keeplimits')
xlabel('UTC','FontSize',14)
ylabel('correlation','FontSize',14)
title(['Hourly Regression Coefficient(' char(name_x) ',' char(name_y) ')'],'FontSize',16)
subplot(1,2,2)
boxplot(r,ceil(mod(utc,1)*24))