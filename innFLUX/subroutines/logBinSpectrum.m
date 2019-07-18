function [f_out, y_out] = logBinSpectrum(f, y, N, f_min, f_max)

bounds(2:(N+1)) = exp( linspace( log(f_min), log(f_max), N ) );
bounds(1) = 0;
bounds(N+1) = f_max;

f_out = NaN(1, N);
y_out = NaN(1, N);

f_out(1) = 0; % first bin contains only the f=0 value
y_out(1) = y(1);
for i = 2:N
    idx = f > bounds(i) & f <= bounds(i+1);
    f_out(i) = exp( (log(bounds(i)) + log(bounds(i+1)) )/2 );
    y_out(i) = mean(y(idx));
end

end
