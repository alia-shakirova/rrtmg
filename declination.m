function delta = declination(time)

% This script calculates declination angle delta in radians
% time is netcdf time
% d is a day number (1 - Jan. 1, 365 - Dec. 31 (366 for leap year))

%time = datetime(double(time)/24 + datenum('1900-01-01 00:00:00'),'ConvertFrom','datenum');
time = double(time)/24 + datenum('1900-01-01 00:00:00');
d = time - datenum(year(time),1,0);
h = hour(time);
%d = day(time,'dayofyear');
y = year(time(1));

if mod(y,4) ~= 0
    theta_d = 2 * pi / 365 * (d - 1 + (h - 12) / 24);
else
    theta_d = 2 * pi / 366 * (d - 1 + (h - 12) / 24);
end
a = [0.006918; -0.399912; -0.006758; -0.002697];
b = [0; 0.070257; 0.000907; 0.001480];
delta = 0;
for n = 0:3
    delta =  delta + a(n+1) * cos(n*theta_d) + b(n+1) * sin(n*theta_d);
end

end
