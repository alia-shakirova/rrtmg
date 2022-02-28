function aver_sza = averaged_zenith_angle(lat,lon,time,step,lim_bot,lim_top)

% This script calculates the averaged cosine solar zenith angle sza 
% time is netcdf time
% lat is latitude (-90, 90)
% lon is longitude
% step is time step in hours
% lim_bot is the bottom limit, lim_top is the top limit

n = (lim_top-lim_bot)/step;
aver_sza = zeros(length(lat),length(lon));

for i = 0:n
    time_new = time - (lim_top - lim_bot) + i * step;
    aver_sza = aver_sza + zenith_angle(lat,lon,time_new);
end

aver_sza = aver_sza/(n+1);
aver_sza(aver_sza>90) = 90;
end
