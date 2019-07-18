function results = innFLUX_qclass(res_in)

%% Classify qaqc flags
%   following Foken et al (2012), Chapter 4 in: Aubinet et al (Eds), Eddy Covariance
%   (see also Foken (2017) DOI 10.1007/978-3-642-25440-6
%    or Foken et al (2004), in: Lee et al (Eds), Handbook of micrometeorology

% (c) martin.graus@uibk.ac.at
% date: 2019-01-08

% Uses  (1) Stationarity test of <w'c'> [SST]
%       (2) Integral turbulence characteristic of u and w [ITC]
%       (3) Flow direction relative to CSAT3 orientation [FD]

% Classes (also adds cmap)
% 1-3 fundamental research (model, spectra,...) [greens]
% 4-6 general use                               [blues]
% 7-8 for orientation only                      [oranges]
% 9 discard for further analysis                [red]
% 10 input was NaN                              [white]

% Combines: (1) into ...qaqc.class(1).val
%           (1)&(2) into ...qaqc.class(2).val
%           (1)&(2)&(3) into ...qaqc.class(3).val
%
% creates statistics of the occurance classes (hist) in each species
% ...qaqc.chist(i).val(1,1:9) to be plotted as stacked bar graph

% assumes res_in being the output of innFLUX_step1 (or _step2)
% and will add fields to res_in (qaqc fields) so can be called by results = innFLUX_qclass(results);
% to just add classification info into original innFLUX results
results = res_in; % write into the original results

% ---MET---
if isfield(results,'MET')
    results.MET.qaqc.SST_wT_class = SST_class(results.MET.qaqc.SST_wT);
    results.MET.qaqc.FlowDir_class = FlowDirection_class(results.MET.wdir,results.parameters.SONIC_ORIENTATION);
    results.MET.qaqc.ITC_uw_class = ITC_class(max(results.MET.qaqc.ITC_w,results.MET.qaqc.ITC_u));
    %   ...qaqc(1).class = indiviual SST
    results.MET.qaqc.class(1).val = results.MET.qaqc.SST_wT_class;  % redundance to keep same sub-STRUCT as MET and TRACER
    %   ...qaqc(2).class = indiviual SST combined with ITC_uw_class
    results.MET.qaqc.class(2).val = ...
        combi_class(results.MET.qaqc.class(1).val,...
        results.MET.qaqc.ITC_uw_class,...
        ones(size(results.MET.qaqc.ITC_uw_class)));
    %   ...qaqc(3).class = indiviual SST combined with ITC_uw_class and FlowDir_class
    results.MET.qaqc.class(3).val = ...
        combi_class(results.MET.qaqc.class(1).val,...
        results.MET.qaqc.ITC_uw_class,...
        results.MET.qaqc.FlowDir_class);
    
    results.MET.qaqc.chist(1).val(1,1:10) = hist(results.MET.qaqc.class(1).val,1:10);
    results.MET.qaqc.chist(2).val(1,1:10) = hist(results.MET.qaqc.class(2).val,1:10);
    results.MET.qaqc.chist(3).val(1,1:10) = hist(results.MET.qaqc.class(3).val,1:10);

    % percentage of class 1-6 in class 1-9 (all but NaN)
    for co_dex = 1:3
    results.MET.qaqc.chist(co_dex).perc_1to6(:,1) = ...
        sum(results.MET.qaqc.chist(co_dex).val(:,1:6),2) ./ ...
        sum(results.MET.qaqc.chist(co_dex).val(:,1:9),2);
    end  
end
% ---end MET---

% ---IRGA---
if isfield(results,'IRGA')
    for I_dex = 1:2 % 1...CO2; 2...H2O
        %   ...qaqc(1).class = indiviual SST
        results.IRGA(I_dex).qaqc.class(1).val = SST_class(results.IRGA(I_dex).qaqc.SST);
        %   ...qaqc(2).class = indiviual SST combined with ITC_uw_class
        results.IRGA(I_dex).qaqc.class(2).val = ...
            combi_class(results.IRGA(I_dex).qaqc.class(1).val,...
            results.MET.qaqc.ITC_uw_class,...
            ones(size(results.MET.qaqc.ITC_uw_class)));
        %   ...qaqc(3).class = indiviual SST combined with ITC_uw_class and FlowDir_class
        results.IRGA(I_dex).qaqc.class(3).val = ...
            combi_class(results.IRGA(I_dex).qaqc.class(1).val,...
            results.MET.qaqc.ITC_uw_class,...
            results.MET.qaqc.FlowDir_class);
        
        results.IRGA(1).qaqc.chist(1).val(I_dex,1:10) = hist(results.IRGA(I_dex).qaqc.class(1).val,1:10);
        results.IRGA(1).qaqc.chist(2).val(I_dex,1:10) = hist(results.IRGA(I_dex).qaqc.class(2).val,1:10);
        results.IRGA(1).qaqc.chist(3).val(I_dex,1:10) = hist(results.IRGA(I_dex).qaqc.class(3).val,1:10);
    end
    
    % percentage of class 1-6 in class 1-9 (all but NaN)
    for co_dex = 1:3
    results.IRGA(1).qaqc.chist(co_dex).perc_1to6(:,1) = ...
        sum(results.IRGA(1).qaqc.chist(co_dex).val(:,1:6),2) ./ ...
        sum(results.IRGA(1).qaqc.chist(co_dex).val(:,1:9),2);
    end    
end
% ---end IRGA---

% VOC (individual and combined with ITC and FD)
if isfield(results,'TRACER')
    for T_dex = 1:length(results.TRACER)
        % individual TRACER
        results.TRACER(T_dex).qaqc.class(1).val = SST_class(results.TRACER(T_dex).qaqc.SST);
        % combined by max(SST_wT, SST_tracer)
        %results.TRACER(T_dex).qaqc.coSST_class = max(results.TRACER(t_dex).qaqc.SST,results.MET.qaqc.SST_wT_class);
        % combined with ITC(u,w) and ...
        results.TRACER(T_dex).qaqc.class(2).val = ...
            combi_class(results.TRACER(T_dex).qaqc.class(1).val,...
            results.MET.qaqc.ITC_uw_class,...
            ones(size(results.MET.qaqc.ITC_uw_class)));
        %... FlowDir_class
        results.TRACER(T_dex).qaqc.class(3).val = ...
            combi_class(results.TRACER(T_dex).qaqc.class(1).val,...
            results.MET.qaqc.ITC_uw_class,...
            results.MET.qaqc.FlowDir_class);
        
        % histogram data of classes
        results.TRACER(1).qaqc.chist(1).val(T_dex,1:10) = hist(results.TRACER(T_dex).qaqc.class(1).val,1:10);
        results.TRACER(1).qaqc.chist(2).val(T_dex,1:10) = hist(results.TRACER(T_dex).qaqc.class(2).val,1:10);
        results.TRACER(1).qaqc.chist(3).val(T_dex,1:10) = hist(results.TRACER(T_dex).qaqc.class(3).val,1:10);
    end

    % percentage of class 1-6 in class 1-9 (all but NaN)
    for co_dex = 1:3
    results.TRACER(1).qaqc.chist(co_dex).perc_1to6(:,1) = ...
        sum(results.TRACER(1).qaqc.chist(co_dex).val(:,1:6),2) ./ ...
        sum(results.TRACER(1).qaqc.chist(co_dex).val(:,1:9),2);
    end
end
% ---end VOC---

% add cmap
results.cmap = [0 0.95 0.3; ...class 1
    0.4 1 0.4;              ...class 2
    0.5 1 0.6;              ...class 3
    0.7 0.85 0.95;          ...class 4
    0.5 0.7 0.95;           ...class 5
    0.3 0.6 1;              ...class 6
    1 0.8 0.3;              ...class 7
    1 0.5 0.2;              ...class 8
    1 0 0;                  ...class 9
    0.95 0.95 0.95];        ...class10 [NaN]
    
