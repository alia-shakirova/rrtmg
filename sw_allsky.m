clear;clc;
tic

  addpath('/home/aliia/RRTM/SW/tests/data');
% parameters
  tind = 1;
  co2 = 403.3; % ppmv % global average co2 in 2016
  addpath('/home/aliia/RRTM/SW/programs/')  
  filename = 'era5_daily_062016_q.nc';
  singfile = 'era5_daily_062016_single.nc';
  outfile    = ['/home/aliia/RRTM/SW/tests/data/hourly/a_0.7_',num2str(tind),'.nc'];
  %outfile    = ['/home/aliia/RRTM/SW/results/ref_',num2str(tind),'.nc'];
  %outfile    = ['/home/aliia/RRTM/SW/results/q_1.1',num2str(tind),'.nc'];
  q_flag = 0;

% read atmospheric data
  lon_coor = ncread(filename,'longitude',1,Inf,4)    % longitude coordinate
  lat_coor = ncread(filename,'latitude',1,Inf,4)    % latitude coordinate
  pmid_tmp = ncread(filename,'level');
  tmid_tmp_in = ncread('era5_daily_062016_t.nc','t',[1 1 1 1],[Inf Inf Inf Inf],[4 4 1 1]);
  h2o_raw_in  = ncread(filename,'q',[1 1 1 1],[Inf Inf Inf Inf],[4 4 1 1]);
  o3_raw_in   = ncread('era5_daily_062016_o3.nc','o3',[1 1 1 1],[Inf Inf Inf Inf],[4 4 1 1]);
  %psfc_raw_in = ncread('sp.nc','sp',[1 1 1],[Inf 61 Inf],[4 4 1]) ./ 100. ;    % surface pressure  % convert to hPa
  psfc_raw_in = ncread(singfile,'sp',[1 1 1],[Inf 61 Inf],[4 4 3]) ./ 100. ; 
  tsfc_raw_in = ncread(singfile,'skt',[1 1 1],[Inf 61 Inf],[4 4 3]);    % skin temperature
  %tisr_raw_in = ncread('tisr.nc','tisr',[1 1 1],[Inf 61 Inf],[4 4 1])./3600;  %TOA incident solar radiation
  tisr_raw_in = ncread(singfile,'tisr',[1 1 1],[Inf 61 Inf],[4 4 3])./3600;
  % add cloud-related varibles
  cloudfrc_in = ncread('era5_daily_062016_cc.nc','cc',[1 1 1 1],[Inf Inf Inf Inf],[4 4 1 1]);    % cloud fraction
  ciwc_raw_in = ncread('era5_daily_062016_ciwc.nc','ciwc',[1 1 1 1],[Inf Inf Inf Inf],[4 4 1 1]);
  clwc_raw_in = ncread('era5_daily_062016_clwc.nc','clwc',[1 1 1 1],[Inf Inf Inf Inf],[4 4 1 1]);
  albedo_raw_in = ncread(singfile,'fal',[1 1 1],[Inf 61 Inf],[4 4 3]);
  tind
  toc
  tic
  tmid_tmp = tmid_tmp_in(:,:,:,tind);
  h2o_raw = h2o_raw_in(:,:,:,tind);
  o3_raw = o3_raw_in(:,:,:,tind);

  psfc_raw = psfc_raw_in(:,:,tind);
  tsfc_raw = tsfc_raw_in(:,:,tind);
  tisr_raw = tisr_raw_in(:,:,tind);
  S0 = max(max(tisr_raw)); %solar constant
  albedo_raw = ones(size(tisr_raw))*0.7;
  %albedo_raw = albedo_raw_in(:,:,tind);
  h2o_raw(h2o_raw < 0) = NaN;        % avoid negative h2o

  cloudfrc = cloudfrc_in(:,:,:,tind);
  ciwc_raw = ciwc_raw_in(:,:,:,tind);
  clwc_raw = clwc_raw_in(:,:,:,tind);
  dim = [length(1:4:1440),length(1:4:121),length(pmid_tmp)+1];

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

% for loop for each gridbox.

  for i =1:dim(1)
  lon = lon_coor(i);
    
    for j = 1:dim(2)
        

        lat = lat_coor(j);
  
        albedo = albedo_raw(i,j);
      
      % calculate hour angle, solar constant
        solar_zenith=acosd(tisr_raw(i,j)/S0);   % retrieved solar zenith angle
 %      solar_zenith = zenith_angle(lat,lon,tind);
       
        if solar_zenith < 90
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
       		 [fup_sw, fdn_dif_sw, fdn_dir_sw, fdn_sw, fnt_sw, htr_sw] = rrtmg_sw_htr_allsky(cflag, atmprofile, solar_zenith, S0, albedo, nlev);
		%size(fup_sw)
		%size(htr_sw)
            else
		fup_sw = 0; fdn_dif_sw = 0; fdn_dir_sw = 0; fdn_sw = 0; fnt_sw = 0; htr_sw = zeros(1,nlev-1);
            end
        	 output_fup_sw(i,j,1:nlev)      = fup_sw;
        	output_fdn_dif_sw(i,j,1:nlev)  = fdn_dif_sw;
        	output_fdn_dir_sw(i,j,1:nlev)  = fdn_dir_sw;
        	output_fdn_sw(i,j,1:nlev)      = fdn_sw;
        	output_fnt_sw(i,j,1:nlev)      = fnt_sw;
        	output_htr_sw(i,j,1:nlev)     = [0 htr_sw];
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
  %dimidt = netcdf.defDim(ncid,'time',tind);
 
  %----------------------------DEFINE NEW VARIABLES-------------------------
 % varid1 = netcdf.defVar(ncid,'fup_sw','NC_FLOAT',[dimidx dimidy dimidz dimidt]);
 % varid2 = netcdf.defVar(ncid,'fdn_dif_sw','NC_FLOAT',[dimidx dimidy dimidz dimidt]);
 % varid3 = netcdf.defVar(ncid,'fdn_dir_sw','NC_FLOAT',[dimidx dimidy dimidz dimidt]);
 % varid4 = netcdf.defVar(ncid,'fdn_sw','NC_FLOAT',[dimidx dimidy dimidz dimidt]);
 % varid5 = netcdf.defVar(ncid,'fnt_sw','NC_FLOAT',[dimidx dimidy dimidz dimidt]);
 % varid6 = netcdf.defVar(ncid,'htr_sw','NC_FLOAT',[dimidx dimidy dimidz dimidt]);
  
  varid1 = netcdf.defVar(ncid,'fup_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid2 = netcdf.defVar(ncid,'fdn_dif_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid3 = netcdf.defVar(ncid,'fdn_dir_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid4 = netcdf.defVar(ncid,'fdn_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid5 = netcdf.defVar(ncid,'fnt_sw','NC_FLOAT',[dimidx dimidy dimidz]);
  varid6 = netcdf.defVar(ncid,'htr_sw','NC_FLOAT',[dimidx dimidy dimidz]);


  %---------------------------DEFINE ATTRIBUTES OF NEW VARIABLES------------
  %netcdf.putAtt(ncid,varid4,'units','hours since 1800-01-01 00:00:0.0');
  
  %---------------------------GIVE VALUES TO VARIABLES-----------------------
  netcdf.endDef(ncid)
  
  netcdf.putVar(ncid,varid1,output_fup_sw);
  netcdf.putVar(ncid,varid2,output_fdn_dif_sw);
  netcdf.putVar(ncid,varid3,output_fdn_dir_sw);
  netcdf.putVar(ncid,varid4,output_fdn_sw);
  netcdf.putVar(ncid,varid5,output_fnt_sw);
  netcdf.putVar(ncid,varid6,output_htr_sw);
  
  netcdf.close(ncid);
  toc
  exit;

