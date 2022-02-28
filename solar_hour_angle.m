function hour_angle = solar_hour_angle(lon,time)

% This script calculates the solar hour angle in radians 
% time is netcdf time
% lon is longitude (positive to the east of the Prime Meridian)
% d is a day number (1 - Jan. 1, 365 - Dec. 31 (366 for leap year))
% Reference: https://www.esrl.noaa.gov/gmd/grad/solcalc/solareqns.PDF

%time = datetime(double(time)/24 + datenum('1900-01-01 00:00:00'),'ConvertFrom','datenum');
time = double(time)/24 + datenum('1900-01-01 00:00:00');
d = time - datenum(year(time),1,0);
sc = second(time);
mn = minute(time);
h = hour(time);
%d = day(time,'dayofyear');
y = year(time(1));

% The fractional year calculation
if mod(y,4) ~= 0
    theta_d = 2 * pi / 365 * (d - 1 + (h - 12) / 24);
else
    theta_d = 2 * pi / 366 * (d - 1 + (h - 12) / 24);  % leap year
end

a = [0.000075; 0.001868; -0.014615];
b = [0; -0.032077; -0.040849];

% Equation of time (in minutes)
eqtime = 0;
for n = 0:2
    eqtime =  eqtime + 229.18 * (a(n+1) * cos(n*theta_d) + b(n+1) * sin(n*theta_d));
end

time_offset = eqtime + 4 * lon;% - 60 * timezone(lon);

% True solar time
tst = h * 60 + mn + sc / 60 + time_offset;
hour_angle = tst / 4 - 180; 
hour_angle = hour_angle * pi / 180; % convert to radians

end
