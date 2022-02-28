function wbroadl = broad(atmosphere,nlev)
    % Calculate column density molecules/cm**2
    % wbroadl 1 x nlev-1
    % from bottom to top
    
    h2o = atmosphere.h2o;
    co2 = atmosphere.co2;
    o3 = atmosphere.o3;
    n2o = atmosphere.n2o;
    co = atmosphere.co;
    ch4 = atmosphere.ch4;
    o2 = atmosphere.o2;
    pint = atmosphere.pint;
    dp = pint(1:nlev-1) - pint(2:nlev);
    
    Na = 6.02e23;
    gravity = 9.8;
    wbroadl = zeros(1,nlev-1);
    
    for ilayer=1:nlev-1
        amm = (1 - h2o(ilayer)) * 28.966 + h2o(ilayer) * 18.016; % The molecular weight of moist air g/mol
        dry_air = dp(ilayer) * 1e3 * Na /(100 * gravity * amm * (1+h2o(ilayer)));
        summol = co2(ilayer) + o3(ilayer) + n2o(ilayer) + co(ilayer) + ch4(ilayer) + o2(ilayer);
        wbroadl(ilayer) = dry_air * (1-summol);
    end