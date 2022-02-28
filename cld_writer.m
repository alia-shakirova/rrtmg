function cld_writer(pint,nlev,cldfrc_tmp,ciwc,clwc)
 
    % make cloud input for RRTNG_SW
        cfile = 'IN_CLD_RRTM';
        fileID=fopen(cfile,'wt');

        inflag = 2;  %0 direct specification of cloud optical depth; 2 calculation separate ice and liquid cloud optical depths
        iceflag = 2; %optical depth
        liqflag = 1; %optical depth

        dp      = pint(1:nlev - 1) - pint(2:nlev);
        indyc   = find(cldfrc_tmp > 0);
        lay     = indyc;
        ncld    = length(indyc);
        cldfrc  = cldfrc_tmp(indyc);
        fracice = ciwc(indyc)./(ciwc(indyc)+clwc(indyc));
    %    effradliq = cldfrc(indyc);
    %    effradice = cldfrc(indyc);
        effradice = 20;
        effradliq = 5;
     
        g = 9.8;

        cwp = 1e2 * 1e3 * (ciwc(indyc) + clwc(indyc)) .* dp(indyc) ./ g ./ (1 + ciwc(indyc) + clwc(indyc));
     
     
        %record C1.1
        fprintf(fileID, '%5i%5i%5i\n', inflag, iceflag, liqflag);
     
        %record C1.2
        for i=1:ncld
            fprintf(fileID, '%5i%10.5f%10.5f%10.5f%10.5f%10.5f\n', lay(i), cldfrc(i), cwp(i), fracice(i), effradice, effradliq);
        end
     
        testchar = '%';
        fprintf(fileID,'%s',testchar);
     
        fclose(fileID);
     
    end
