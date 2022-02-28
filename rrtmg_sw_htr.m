function [fup, fdn_dif, fdn_dir, fdn, fnt, htr] = rrtmg_sw_htr(atmprofile, solar_zenith, S0, albedo, nlev)

    ! rm -rf TAPE5,TAPE6,OUTPUT_RRTM
    
      rrtmg_tape5_writer_htr_sw(atmprofile, solar_zenith, S0, albedo, nlev);
    
    !./rrtmg_sw
    
      [fup, fdn_dif, fdn_dir, fdn, fnt, htr1] = rrtmg_sw_output_read('OUTPUT_RRTM',nlev);
    htr = htr1(2:nlev);
