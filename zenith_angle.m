function sza = zenith_angle(lat,lon,time)

% This script calculates the cosine solar zenith angle sza 
% time is netcdf time
% delta is the declination angle
% h is hour angle
% lat is latitude (-90, 90)
% lon is longitude

if lon > 180
    lon = lon - 360;
end
lat = lat * pi / 180;
h = solar_hour_angle(lon,time);
delta = declination(time);

sza = acos(sin(lat) * sin(delta) * ones(size(h))' + cos(lat) * cos(delta) * cos(h)');
sza = sza / pi * 180;
sza(sza>90) = 90;
end


