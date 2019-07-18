function B_classDOW = innFLUX_DayOfWeek_class(holidays,utc)
%% Classification of Day of Week (DOW) by boolean variable accounting for holidays
% add meta to STRUCT
b(1).meta ={...
    'b(d_dex).day = 1...7 (Sun, Mon-Sat AND NOT[Holiday])'; ...
	'b(1).Hol = 1 <=> Holiday'; ...
	'b(1).SunHol = 1 <=> Sunday OR Holiday'; ...
	'b(1).MonFri = 1 <=> Mon-Fri AND NOT[Holiday]'; ...
	'b(1).TueThu = 1 <=> Tue-Thu AND NOT[Holiday]'};

% boolian indicating holidays (local time)
for h_dex = 1:7
    ddex_h(1:length(utc),h_dex) = (floor(utc) == holidays(h_dex));
end
b(1).hol = sum(ddex_h,2,'native');

% sunday 
d_dex = 1;
b(d_dex).day = (weekday(utc)==d_dex);
% workday AND NOT holiday
for d_dex = 2:7 
    b(d_dex).day = (~(b(1).hol) & weekday(utc)==d_dex);
end

% any workday (Mon-Fri) AND NOT holiday
b(1).MonFri = sum(cat(2,b(2).day,b(3).day,b(4).day,b(5).day,b(6).day),2);

% Tue-Thu AND NOT holiday
b(1).TueThu = sum(cat(2,b(3).day,b(4).day,b(5).day),2);

% copy of Sun OR holiday
b(1).SunHol = (b(1).day | b(1).hol); 

% boolean by month
for m_dex = 1:12
    b(1).month(1:length(utc),m_dex) = (month(utc) == m_dex);
end

B_classDOW = b;
return

end
