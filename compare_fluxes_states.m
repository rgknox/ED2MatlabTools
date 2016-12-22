%==========================================================================
% Program: compare_fluxes_states.m
%
% This script runs a suit of checks on time-series ED2 output.
% The output variables are from "Q" files, and also are expected to be
% polygon scale.  Do not try to use other variables.
% The user can provide data streams for up to four different simulation
% outputs, defined by "dset".
%
%==========================================================================

clear all;
close all;

%==========================================================================
%     User defined variables
%==========================================================================

test_name = 'ZF2';
figpref   = 'zf2_v1';

dset(1).lab = '1PCNT';
dset(1).Spfx = '';
dset(1).Qpfx = 'TestData/ZF2_TF1PCNT/zf2_oxi_TF1PCNT-Q';

dset(2).lab = '.9PCNT';
dset(2).Spfx = '';
dset(2).Qpfx = 'TestData/ZF2_TF.9PCNT/zf2_oxi_TF.9PCNT-Q';


fluxvar(1).printname = 'ET';
fluxvar(1).varname = '/QMEAN_VAPOR_AC_PY';
fluxvar(1).units = '[mm/m^2/mo]';
fluxvar(1).scalar = -86400*30;
fluxvar(1).offset = 0.0;
fluxvar(1).is2d = false;

fluxvar(2).printname = 'SHF';
fluxvar(2).varname = '/QMEAN_SENSIBLE_AC_PY';
fluxvar(2).units = '[W/m^2]';
fluxvar(2).scalar = 1.0;
fluxvar(2).offset = 0.0;
fluxvar(2).is2d = false;

fluxvar(3).printname = 'Rnet';
fluxvar(3).varname = '/QMEAN_RNET_PY';
fluxvar(3).units = '[W/m^2]';
fluxvar(3).scalar = 1.0;
fluxvar(3).offset = 0.0;
fluxvar(3).is2d = false;

fluxvar(4).printname = 'Rswd';
fluxvar(4).varname = '/QMEAN_RSHORT_L_PY';
fluxvar(4).units = '[W/m^2]';
fluxvar(4).scalar = 1.0;
fluxvar(4).offset = 0.0;
fluxvar(4).is2d = false;

fluxvar(5).printname = 'Rlwd';
fluxvar(5).varname = '/QMEAN_RLONG_L_PY';
fluxvar(5).units = '[W/m^2]';
fluxvar(5).scalar = 1.0;
fluxvar(5).offset = 0.0;
fluxvar(5).is2d = false;

fluxvar(6).printname = 'GPP';
fluxvar(6).varname = '/QMEAN_GPP_PY';
fluxvar(6).units = '[kgC/m^2/yr]';
fluxvar(6).scalar = 1.0;
fluxvar(6).offset = 0.0;
fluxvar(6).is2d = false;

fluxvar(7).printname = 'NEP';
fluxvar(7).varname = '/QMEAN_NEP_PY';
fluxvar(7).units = '[kgC/m^2/yr]';
fluxvar(7).scalar = 1.0;
fluxvar(7).offset = 0.0;
fluxvar(7).is2d = false;

fluxvar(8).printname = 'Canopy CO2';
fluxvar(8).varname = '/QMEAN_CAN_CO2_PY';
fluxvar(8).units = '[ppm]';
fluxvar(8).scalar = 1.0;
fluxvar(8).offset = 0.0;
fluxvar(8).is2d = false;

fluxvar(9).printname = 'Mean Soil Moisture';
fluxvar(9).varname = '/QMEAN_SOIL_WATER_PY';
fluxvar(9).units = '[m^3/m^3]';
fluxvar(9).scalar = 1.0;
fluxvar(9).offset = 0.0;
fluxvar(9).is2d = true;

fluxvar(10).printname = 'Leaf Temp';
fluxvar(10).varname = '/QMEAN_LEAF_TEMP_PY';
fluxvar(10).units = '[^oC]';
fluxvar(10).scalar = 1.0;
fluxvar(10).offset = -273.14;
fluxvar(10).is2d = false;

nfvar = length(fluxvar(:));

if(length(dset(:))>4)
    display('Can only intercompare up to 4 simulations');
    return
end

if(nfvar==0)
    display('Check your file prefix.');
    display('We didnt find any files with that prefix.');
    display('Note also that the word didnt in that last sentence');
    display('should have a single apostrophe.  But... adding that apostrophe')
    display('would had been overly difficult given that single apostrophes')
    display('are used to dilineate text strings in matlab');
    display('So please just imagine that words that seem like they should');
    display('have appostrophes... do have apostrophes.');
    return
end


global fasz;
fasz = 12;
visible = 'on';
outdir = sprintf('Figures/',test_name);
addpath('Functions/');

ndsets = numel(dset);


% First, force the sample sizes to be the same
% ============================================

minqfiles = 1e10;
for kd=1:ndsets
    Qlist = dir(strcat(dset(kd).Qpfx,'*h5'));
    id   = strfind(dset(kd).Qpfx,'/');
    Qdir = dset(kd).Qpfx(1:id(end));
    nqfiles     = length(Qlist);
    minyear = 1e10;
    maxyear = -1e10;
    for it=1:nqfiles
        qfile = strcat(Qdir,Qlist(it).name);
        iyear  = str2double(qfile(end-23:end-20));
        maxyear = max(iyear,maxyear);
        minyear = min(iyear,minyear);
    end
    dset(kd).minyear=minyear;
    dset(kd).maxyear=maxyear;
    minqfiles = min([nqfiles,minqfiles]);
end

for kd=1:ndsets
    
    Qlist = dir(strcat(dset(kd).Qpfx,'*h5'));
    Qlist = Qlist(end-minqfiles+1:end);
    
    id   = strfind(dset(kd).Qpfx,'/');
    Qdir = dset(kd).Qpfx(1:id(end));
    
    nqfiles     = length(Qlist);
    dnq = zeros(nqfiles,1);
    
    %==================================================================
    % Comparison of Flux Variables
    %==================================================================
    
    for it=1:nqfiles
        
        qfile = strcat(Qdir,Qlist(it).name);
        
        iyear  = str2double(qfile(end-23:end-20));
        imonth = str2double(qfile(end-18:end-17));
        idate  = str2double(qfile(end-15:end-14));
        ihour  = str2double(qfile(end-12:end-11));
        iminute= str2double(qfile(end-10:end-9));
        isecond= str2double(qfile(end-8:end-7));
        
        dnq(it) = datenum(iyear,imonth,idate,ihour,iminute,isecond);
        
        if(it==nqfiles)
            iyearz=iyear;
        end
        
        if(it==1)
            iyeara=iyear;
            nhrs = numel(hdf5read(qfile,'/QMEAN_VAPOR_AC_PY'));
            tmp = -hdf5read(qfile,'/SLZ');
            k50 = find(tmp>0.5,1,'last');
            nz = numel(tmp);
            dz = zeros(nz,1);
            for iz=1:nz-1
                dz(iz)=tmp(iz)-tmp(iz+1);
            end
            dz(nz) = tmp(nz);
            slz=zeros(nz,nhrs);
            for ihr=1:nhrs
                slz(:,ihr)=dz;
            end
            
            phrs = linspace(24/nhrs,24,nhrs);
            npts = nqfiles*nhrs;
            nmos    = zeros(12,1);
            
            flux_tp = zeros(npts,nfvar);
            flux_td = zeros(nhrs,nfvar);
            flux_tm = zeros(12,nfvar);
            flux_cp = zeros(npts,nfvar);
            flux_cd = zeros(nhrs,nfvar);
            flux_cm = zeros(12,nfvar);
            
        end
        
        id1 = (it-1)*nhrs+1;
        id2 = (it*nhrs);
        nmos(imonth)=nmos(imonth)+1;
        
        for itmp = 1:nfvar
            if(fluxvar(itmp).is2d)
                tmp2d = hdf5read(qfile,fluxvar(itmp).varname).*fluxvar(itmp).scalar + fluxvar(itmp).offset;
                tmp   = sum(tmp2d(:,:).*slz(:,:),1)'./ sum(slz(:,1));
            else
                tmp = hdf5read(qfile,fluxvar(itmp).varname)*fluxvar(itmp).scalar + fluxvar(itmp).offset;
            end
            flux_tp(id1:id2,itmp) = tmp;
            flux_td(:,itmp)      = flux_td(:,itmp) + tmp./nqfiles;
            flux_tm(imonth,itmp) = flux_tm(imonth,itmp)+mean(tmp);
            flux_names{itmp} = fluxvar(itmp).printname;
            flux_units{itmp} = fluxvar(itmp).units;
        end
        
    end
    
    % Normalize means
    for imo=1:12
        flux_tm(imo,:)  = flux_tm(imo,:)./nmos(imo);
    end
    
    dset(kd).flux_tm = flux_tm;
    dset(kd).flux_td = flux_td;
    dset(kd).flux_tp = flux_tp;
   
end


% Make a plot
%==================================================================
fluxes_img = sprintf('%sfluxes_%s.png',outdir,figpref);
titlestr = sprintf('%s (%4i - %4i)',test_name,iyeara,iyearz);
plot_fluxes(dset,flux_names,flux_units,phrs,fluxes_img,visible)





