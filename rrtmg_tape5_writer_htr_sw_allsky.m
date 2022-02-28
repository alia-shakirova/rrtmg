
function [] = rrtmg_tape5_writer_htr_sw_allsky(cflag,atmprofile,solar_zenith, S0, albedo, nlev)

    filename = 'INPUT_RRTM';
    fileID=fopen(filename,'wt');
    
    %shortwave set record 1.1~basic setting
    record1.iaer= 0;		%flag for aerosols (0,10) 0 no  aerosols 6 ECMWF global mean aerosol properties 10 user defined IN_AER_RRTM
    record1.iatm     = 0;		%flag for RRTATM (0,1->yes)
    record1.iscat    = 1;
    record1.istrm    = 0;
    record1.iout     = 0;		%-1 no output% 0 for 10-3250% n from band n% 99 17 spectral bands
    record1.imca     = cflag.imca;		%flag for McICA of sub-grid cloud fraction
    record1.icld     = cflag.icld;		%0-no cloudy%1-RANDOM overlap assmption(IMAC=0 or 1)%2-MAXIMUM/RANDUM(IMAC=0 or 1)% 3-MAXIMUM(IMAC=1)
    record1.idelm    = 0;
    record1.icos     = 0;
    
    %solar setting
    record1.juldat   = 0;		%0-no Earth-sun distance variation
    record1.sza      = solar_zenith;
    record1.isolvar  = 0;
    record1.scon     = S0;
    
    %surface setting
    record1.iemis    = 1;		%surface emissivity (0-1.0,1-same emissivity,2-different)
    record1.ireflect = 0;	
    record1.semiss   = 1 - albedo;	%surface emissivity
    
    %atmospheric profile
    record2.iform    = 0;
    record2.nlayrs   = nlev - 1;
    record2.nmol     = 7;
    record2.pave     = atmprofile.pmid;
    record2.tave     = atmprofile.tmid;
    record2.pz       = atmprofile.pint;
    record2.tz       = atmprofile.tint;
    record2.wkl(1,:) = atmprofile.h2o;
    record2.wkl(2,:) = atmprofile.co2;
    record2.wkl(3,:) = atmprofile.o3;
    record2.wkl(4,:) = atmprofile.n2o;
    record2.wkl(5,:) = atmprofile.co;
    record2.wkl(6,:) = atmprofile.ch4;
    record2.wkl(7,:) = atmprofile.o2;
    record2.wbroadl  = atmprofile.wbroadl;
    
    record1.slist    = '$ ATMOSPHERE PROFILE';
  
    %record 1.1
    fprintf(fileID, '%s\n',record1.slist);
  
    %record 1.2
    fprintf(fileID,'%20d%30d%33d%2d%5d%4d%1d%4d%1d\n',record1.iaer, record1.iatm,...
        record1.iscat, record1.istrm, record1.iout, record1.imca,...
        record1.icld, record1.idelm, record1.icos);
    %record 1.2.1
  
    fprintf(fileID,'%15d%10.4f%5d%10.4f\n',record1.juldat, record1.sza,...
        record1.isolvar, record1.scon);
    %record 1.4
  
    fprintf(fileID,'%12d%3d%5.3f\n',record1.iemis, record1.ireflect, record1.semiss);
  
    % record 2.1
    fprintf(fileID,'%2d%3d%5d\n',record2.iform, record2.nlayrs, record2.nmol);
  
    % record 2.1.1
    for ilayer=1:nlev-1
        fprintf(fileID,'%10.4f%10.4f%31.3f%7.2f%15.3f%7.2f\n',record2.pave(ilayer), record2.tave(ilayer),...
            record2.pz(ilayer), record2.tz(ilayer), record2.pz(ilayer+1), record2.tz(ilayer+1));
        fprintf(fileID,'%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e%10.3e\n',record2.wkl(1,ilayer),...
            record2.wkl(2,ilayer), record2.wkl(3,ilayer), record2.wkl(4,ilayer), record2.wkl(5,ilayer),...
            record2.wkl(6,ilayer), record2.wkl(7,ilayer), record2.wbroadl(ilayer));
    end
  
    fprintf(fileID, '%s','%%%%%');
    fclose(fileID);
  
  end
