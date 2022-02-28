clear;clc;
tic

  addpath('/home/aliia/RRTM/SW/dataset/');
% parameters
  tind = 1;
  co2 = 400; % ppmv % global average co2 
  addpath('/home/aliia/RRTM/SW/script/')  
  time_period = '062016';  % month year
  aind = 1.0;
  filename = ['era5_hourly_',time_period,'_pressurelevels.nc'];
  singfile = ['era5_hourly_',time_period,'_single.nc'];
  outfile    = ['/home/aliia/RRTM/SW/results/albedo_feedback/',time_period,'_a',num2str(aind),'_',num2str(tind),'.nc'];
  outfile2    = ['/home/aliia/RRTM/SW/results/albedo_feedback/',time_period,'_a',num2str(aind),'_',num2str(tind),'_toa_sfc.nc'];
  q_flag = 0;

% read atmospheric data
  minlat = 21; % 70N
  lon_coor = ncread(filename,'longitude',1,Inf);    % longitude coordinate
  lat_coor = ncread(filename,'latitude',1,minlat);    % latitude coordinate
  pmid_tmp = ncread(filename,'level');
  tmid_tmp_in = ncread(filename,'t',[1 1 1 1],[Inf minlat Inf Inf],[1 1 1 1]);
  h2o_raw_in  = ncread(filename,'q',[1 1 1 1],[Inf minlat Inf Inf],[1 1 1 1]);
  o3_raw_in   = ncread(filename,'o3',[1 1 1 1],[Inf minlat Inf Inf],[1 1 1 1]);
  time = ncread(filename,'time');
  time = time(tind);
  psfc_raw_in = ncread(singfile,'sp',[1 1 1],[Inf minlat Inf],[1 1 1]) ./ 100. ;
  tsfc_raw_in = ncread(singfile,'skt',[1 1 1],[Inf minlat Inf],[1 1 1]);    % skin temperature
  %tisr_raw_in = ncread(singfile,'tisr',[1 1 1],[Inf minlat Inf],[1 1 1])./3600;  %TOA incident solar radiation
  % add cloud-related varibles
  cloudfrc_in = ncread(filename,'cc',[1 1 1 1],[Inf minlat Inf Inf],[1 1 1 1]);    % cloud fraction
  ciwc_raw_in = ncread(filename,'ciwc',[1 1 1 1],[Inf minlat Inf Inf],[1 1 1 1]);
  clwc_raw_in = ncread(filename,'clwc',[1 1 1 1],[Inf minlat Inf Inf],[1 1 1 1]);
  albedo_raw_in = ncread(singfile,'fal',[1 1 1],[Inf minlat Inf],[1 1 1]);
  tind
  toc
  tic
  tmid_tmp = tmid_tmp_in(:,:,:,tind);
  h2o_raw = h2o_raw_in(:,:,:,tind);
  o3_raw = o3_raw_in(:,:,:,tind);

  psfc_raw = psfc_raw_in(:,:,tind);
  tsfc_raw = tsfc_raw_in(:,:,tind);
  tisr_raw = tisr_raw_in(:,:,tind);
  %albedo_raw = albedo_raw_in(:,:,tind);
  albedo_raw = aind * ones(size(psfc_raw));
  albedo_raw(1,1)
  h2o_raw(h2o_raw < 0) = NaN;        % avoid negative h2o

  time_new = double(time)/24 + datenum('1900-01-01 00:00:00');
  d = time_new - datenum(year(time_new),1,0);
  y = year(time_new);
  h = hour(time_new);

  % The fractional year calculation
  if mod(y,4) ~= 0
      theta_d = 2 * pi / 365 * (d - 1 + (h - 12) / 24);
  else
      theta_d = 2 * pi / 366 * (d - 1 + (h - 12) / 24);  % leap year
  end
  coef = 1.000110+0.034221*cos(theta_d)+0.001280*sin(theta_d)+0.000719*cos(2*theta_d)+0.000077*sin(2*theta_d);  
  S0 = 1360 * coef; % solar constant

  cloudfrc = cloudfrc_in(:,:,:,tind);
  ciwc_raw = ciwc_raw_in(:,:,:,tind);
  clwc_raw = clwc_raw_in(:,:,:,tind);
  dim = [length(1:360),length(1:minlat),length(pmid_tmp)+1];
  clear cloudfrc_in ciwc_raw_in clwc_raw_in tsfc_raw_in tisr_raw_in psfc_raw_in tmid_tmp_in h2o_raw_in o3_raw_in albedo_raw_in
  path_rrtmg_sw = '/home/aliia/rrtmg_sw_v4.0/column_model/build/';
  rrtmg_sw_exe = 'rrtmg_sw_v3.9_linux_intel'; 
  unix(['ln -sf ',path_rrtmg_sw, rrtmg_sw_exe, ' ./rrtmg_sw']);

% Create struct for output data
  output_fup_sw     =  NaN(dim);
  output_fdn_dif_sw = NaN(dim);
  output_fdn_dir_sw = NaN(dim);
  output_fdn_sw    = NaN(dim);
  output_fnt_sw   = NaN(dim);
  output_htr_sw    = NaN(dim);
  fup_sw_toa = NaN(dim(1:2));
  fup_sw_sfc = NaN(dim(1:2));
  fdn_sw_toa = NaN(dim(1:2));
  fdn_sw_sfc = NaN(dim(1:2));
  fdn_dif_sw_toa = NaN(dim(1:2));
  fdn_dif_sw_sfc = NaN(dim(1:2));
  fdn_dir_sw_toa = NaN(dim(1:2));
  fdn_dir_sw_sfc = NaN(dim(1:2));
  fnt_sw_toa = NaN(dim(1:2));
  fnt_sw_sfc = NaN(dim(1:2));
% for loop for each gridbox.

  for i =1:dim(1)
    lon = lon_coor(i);
    
    for j = 1:dim(2)
        
      lat = lat_coor(j);
      albedo = albedo_raw(i,j);
    
      % calculate hour angle, solar constant
      sza = averaged_zenith_angle(lat,lon,time,0.01,-1,0);
     
      if sza < 90
        psfc = psfc_raw(i,j);
        tsfc = tsfc_raw(i,j);
      
        if q_flag == 1 
          if j <= 61
               h2o_raw(i,j,18:end) = h2o_raw(i,j,18:end)*1.1;
          else
              h2o_raw(i,j,17:end) = h2o_raw(i,j,17:end)*1.1;
          end
        end

         % read profile
        [nlev, pmid, pint, tmid, tint, n_to_cut, h2o, o3, n2o, co, ch4, o2, cflag, cldfrc_tmp] = read_profile_allsky(filename,i,j,pmid_tmp,tmid_tmp,psfc,tsfc,h2o_raw,o3_raw,cloudfrc);
        
        atmprofile.pmid    = pmid;
        atmprofile.pint    = pint;
        atmprofile.tmid    = tmid;
        atmprofile.tint    = tint;
        atmprofile.h2o     = h2o;
        atmprofile.o3      = o3;
        atmprofile.n2o     = n2o;
        atmprofile.co      = co;
        atmprofile.ch4     = ch4;
        atmprofile.o2      = o2;
        atmprofile.co2     = co2.*ones(size(atmprofile.h2o)) * 1e-6;
        atmprofile.wbroadl = broad(atmprofile,nlev);

        if cflag.imca == 1
          ciwc = flip(squeeze(ciwc_raw(i,j,:)))';
          clwc = flip(squeeze(clwc_raw(i,j,:)))';
          ciwc   = ciwc(n_to_cut+1:dim(3)-1);
          clwc   = clwc(n_to_cut+1:dim(3)-1);
          cld_writer(pint,nlev,cldfrc_tmp,ciwc,clwc);
        end
      
      %disp('calculate sw')
        % calculate shortwave
        [fup_sw, fdn_dif_sw, fdn_dir_sw, fdn_sw, fnt_sw, htr_sw] = rrtmg_sw_htr_allsky(cflag, atmprofile, sza, S0, albedo, nlev);
      else
        nlev = length(pmid_tmp);
        fup_sw = 0; 
        fdn_dif_sw = 0; 
        fdn_dir_sw = 0; 
        fdn_sw = 0; 
        fnt_sw = 0; 
        htr_sw = zeros(1,nlev-1);
      end
      output_fup_sw(i,j,1:nlev)      = fup_sw;
      output_fdn_dif_sw(i,j,1:nlev)  = fdn_dif_sw;
      output_fdn_dir_sw(i,j,1:nlev)  = fdn_dir_sw;
      output_fdn_sw(i,j,1:nlev)      = fdn_sw;
      output_fnt_sw(i,j,1:nlev)      = fnt_sw;
      output_htr_sw(i,j,1:nlev)     = [0 htr_sw];
      fup_sw_toa(i,j) = fup_sw(1);
      fup_sw_sfc(i,j) = fup_sw(find(~isnan(fup_sw),1,'last'));
      fdn_sw_toa(i,j) = fdn_sw(1);
      fdn_sw_sfc(i,j) = fdn_sw(find(~isnan(fdn_sw),1,'last'));
      fdn_dif_sw_toa(i,j) = fdn_dif_sw(1);
      fdn_dif_sw_sfc(i,j) = fdn_dif_sw(find(~isnan(fdn_dif_sw),1,'last'));
      fdn_dir_sw_toa(i,j) = fdn_dir_sw(1);
      fdn_dir_sw_sfc(i,j) = fdn_dir_sw(find(~isnan(fdn_dir_sw),1,'last'));
      fnt_sw_toa(i,j) = fnt_sw(1);
      fnt_sw_sfc(i,j) = fnt_sw(find(~isnan(fnt_sw),1,'last'));
    end
  end
 toc

 tic
 disp('write to nc file')
  ! rm -rf outfile.nc
  
  % --------------------------- DEFINE THE FILE ----------------------------
  ncid = netcdf.create(outfile,'CLOBBER');
  
  %-----------------------------DEFINE DIMENSION----------------------------
  dimidx = netcdf.defDim(ncid,'lon',dim(1));
  dimidy = netcdf.defDim(ncid,'lat',dim(2));
  dimidz = netcdf.defDim(ncid,'level',dim(3));
 
  %----------------------------DEFINE NEW VARIABLES-------------------------
  varid1 = netcdf.defVar(ncid,'fup_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid2 = netcdf.defVar(ncid,'fdn_dif_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid3 = netcdf.defVar(ncid,'fdn_dir_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid4 = netcdf.defVar(ncid,'fdn_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid5 = netcdf.defVar(ncid,'fnt_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid6 = netcdf.defVar(ncid,'htr_sw','NC_FLOAT',[dimidx dimidy dimidz]);

  %---------------------------GIVE VALUES TO VARIABLES-----------------------
  netcdf.endDef(ncid)
  netcdf.putVar(ncid,varid1,output_fup_sw);
  netcdf.putVar(ncid,varid2,output_fdn_dif_sw);
  netcdf.putVar(ncid,varid3,output_fdn_dir_sw);
  netcdf.putVar(ncid,varid4,output_fdn_sw);
  netcdf.putVar(ncid,varid5,output_fnt_sw);
  netcdf.putVar(ncid,varid6,output_htr_sw);
  
  netcdf.close(ncid);

  disp('write to nc file')
  ! rm -rf outfile2.nc
  
  % --------------------------- DEFINE THE FILE ----------------------------
  ncid2 = netcdf.create(outfile2,'CLOBBER');
  
  %-----------------------------DEFINE DIMENSION----------------------------
  dimidx = netcdf.defDim(ncid2,'lon',dim(1));
  dimidy = netcdf.defDim(ncid2,'lat',dim(2));
 % dimidz = netcdf.defDim(ncid2,'level',dim(3));
 
  %----------------------------DEFINE NEW VARIABLES-------------------------
  varid1 = netcdf.defVar(ncid2,'fup_sw_toa','NC_FLOAT',[dimidx dimidy]);
  varid2 = netcdf.defVar(ncid2,'fup_sw_sfc','NC_FLOAT',[dimidx dimidy]);
  varid3 = netcdf.defVar(ncid2,'fdn_dif_sw_toa','NC_FLOAT',[dimidx dimidy]);
  varid4 = netcdf.defVar(ncid2,'fdn_dif_sw_sfc','NC_FLOAT',[dimidx dimidy]);
  varid5 = netcdf.defVar(ncid2,'fdn_dir_sw_toa','NC_FLOAT',[dimidx dimidy]);
  varid6 = netcdf.defVar(ncid2,'fdn_dir_sw_sfc','NC_FLOAT',[dimidx dimidy]);
  varid7 = netcdf.defVar(ncid2,'fdn_sw_toa','NC_FLOAT',[dimidx dimidy]);
  varid8 = netcdf.defVar(ncid2,'fdn_sw_sfc','NC_FLOAT',[dimidx dimidy]);
  varid9 = netcdf.defVar(ncid2,'fnt_sw_toa','NC_FLOAT',[dimidx dimidy]);
  varid10 = netcdf.defVar(ncid2,'fnt_sw_sfc','NC_FLOAT',[dimidx dimidy]);

  %---------------------------GIVE VALUES TO VARIABLES-----------------------
  netcdf.endDef(ncid2)
  netcdf.putVar(ncid2,varid1,fup_sw_toa);
  netcdf.putVar(ncid2,varid2,fup_sw_sfc);
  netcdf.putVar(ncid2,varid3,fdn_dif_sw_toa);
  netcdf.putVar(ncid2,varid4,fdn_dif_sw_sfc);   
  netcdf.putVar(ncid2,varid5,fdn_dir_sw_toa);
  netcdf.putVar(ncid2,varid6,fdn_dir_sw_sfc);
  netcdf.putVar(ncid2,varid7,fdn_sw_toa);
  netcdf.putVar(ncid2,varid8,fdn_sw_sfc);
  netcdf.putVar(ncid2,varid9,fnt_sw_toa);
  netcdf.putVar(ncid2,varid10,fnt_sw_sfc);
  
  netcdf.close(ncid2);
  toc
  exit;

