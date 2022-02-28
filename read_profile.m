function [nlev, pmid, pint, tmid, tint, h2o, o3, n2o, co, ch4, o2] = read_profile(filename,i,j,pmid_tmp,tmid_tmp,psfc,tsfc,h2o_raw,o3_raw)

    pmid_tmp = cast(pmid_tmp,'single');	% convert int32 to single for interp
    pmid_tmp = flip(pmid_tmp)';	% reverse pmid_tmp to bottom-up % make it consistent with example file		% transpose matrix
    nlev = length(pmid_tmp); 	% will be updated to actual layers
    
    % cut unrealistic layers
    n_to_cut = length(find(pmid_tmp > psfc));
    
    pmid = pmid_tmp(n_to_cut+1:nlev);
    
    tmid_tmp = flip(squeeze(tmid_tmp(i,j,:)))';
    tmid_tmp = cast(tmid_tmp,'single');  %test float
    tmid = tmid_tmp(n_to_cut+1:nlev);
    
    % add surface layer
    pint = [psfc, (pmid(1:length(pmid)-1)+pmid(2:length(pmid)))./2, 0.5];	% set model top at 0.5hPa
    tint = [tsfc, interp1(log(pmid), tmid, log(pint(2:length(pint))), 'linear', 'extrap')];
    
    h2o = flip(squeeze(h2o_raw(i,j,:)))';
    h2o = cast(h2o,'single')* 28.960./ 18.016;  % specific humidity -> volume mixing ratio
    h2o = h2o(n_to_cut+1:nlev);
    
    o3  = flip(squeeze(o3_raw(i,j,:)))';
    o3  = cast(o3,'single')* 28.960 ./47.998; 	% mass mixing ratio -> volume mixing ratio
    o3  = o3(n_to_cut+1:nlev);
    
    % use the latest n2o and ch4 concentration
    % https://public.wmo.int/en/media/press-release/greenhouse-gas-levels-atmosphere-reach-new-record
    n2o = 0.329 * ones(1,length(h2o)) * 1e-6;
    co  = 0.15  * ones(1,length(h2o)) * 1e-6;
    ch4 = 1.859 * ones(1,length(h2o)) * 1e-6;
    o2  = 0.209 * ones(1,length(h2o));
    
    % update nlev after cutting layers
    nlev = length(pint);
    
end
