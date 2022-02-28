function [fup, fdn_dif, fdn_dir, fdn, fnt, htr] = rrtmg_sw_output_read(filename,nlev)

    htr = zeros(1,nlev);
    fnt = zeros(1,nlev);
    fup = zeros(1,nlev);
    fdn_dif = zeros(1,nlev);
    fdn_dir = zeros(1,nlev);
    fdn = zeros(1,nlev);
    fi = fopen(filename,'r');
    for j = 1:5; tmp1 = fgetl(fi); end % read through the header info
    for j =1:nlev
        tmp1 = fgetl(fi);
        tmp2 = sscanf(tmp1,'%g',[8 1]);
        %     p(j) = tmp2(2);
        fup(j) = tmp2(3);
        fdn_dif(j) = tmp2(4);
        fdn_dir(j) = tmp2(5);
        fdn(j) = tmp2(6);
        fnt(j) = tmp2(7);
        htr(j) = tmp2(8);
    end
    fclose(fi);
