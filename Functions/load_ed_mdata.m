function [mdata]=load_ed_mdata(edqpfx)


% =========================================================================
% Load ED Model output, using "Q"-files
% =========================================================================

Qlist = dir(strcat(edqpfx,'*h5'));

id   = strfind(edqpfx,'/');
Qdir = edqpfx(1:id(end));

nqfiles = length(Qlist);
dnq     = zeros(nqfiles,1);


% Find out how many patches we have total, in all q files
% The assumption is that you must have all months for each year
%==========================================================================

iyear   = zeros(nqfiles,1);
imonth  = zeros(nqfiles,1);
npatchq = zeros(nqfiles,1);
areaq   = zeros(nqfiles,1);

for it=1:nqfiles
    
    edfile = strcat(Qdir,Qlist(it).name);
    iyear(it)  = str2double(edfile(end-23:end-20));
    imonth(it) = str2double(edfile(end-18:end-17));
    pysi_id = hdf5read(edfile,'/PYSI_ID');
    pysi_n  = hdf5read(edfile,'/PYSI_N');
    sipa_id = hdf5read(edfile,'/SIPA_ID');
    sipa_n  = hdf5read(edfile,'/SIPA_N');
    dist_type  = hdf5read(edfile,'/DIST_TYPE');
    patch_area = hdf5read(edfile,'/AREA');
    site_area  = hdf5read(edfile,'/AREA_SI');
    isi_b = pysi_id(1);
    isi_e = isi_b+pysi_n(1) - 1;
    
    % Count patches and force area to unity
    npatch=0;
    for isi=isi_b:isi_e
        ipa_b = sipa_id(isi);
        ipa_e = ipa_b + sipa_n(isi) - 1;
        for ipa=ipa_b:ipa_e
            if dist_type(ipa)==3 % Primary Forest Only
                npatch=npatch+1;
                areaq(it) = areaq(it)+site_area(isi).*patch_area(ipa);
            end
        end
    end
    npatchq(it) = npatch;
end

% Sanity Check to make sure you have all months for each year of interst
maxyear = max(iyear);
minyear = min(iyear);
for iy=minyear:maxyear
    ids=find(iyear==iy);
    if(numel(ids)~=12)
        display('YOU MUST HAVE ALL 12 MONTHS PER EACH YEAR OF Q FILES');
        return;
    end
end

nyears_ed = maxyear-minyear+1;

max_dbh_early = 0;
max_dbh_late  = 0;

for it=1:nqfiles
    
    edfile = strcat(Qdir,Qlist(it).name);
    
    iy  = str2double(edfile(end-23:end-20))-minyear+1;
    im  = str2double(edfile(end-18:end-17));
    
    mdata.date(iy,im).min_npaf_all = 100;
    
    mdata.date(iy,im).year = str2double(edfile(end-23:end-20));
    
    dbhmin = 0.0;
    
    pysi_id = hdf5read(edfile,'/PYSI_ID');
    pysi_n  = hdf5read(edfile,'/PYSI_N');
    sipa_id = hdf5read(edfile,'/SIPA_ID');
    sipa_n  = hdf5read(edfile,'/SIPA_N');
    paco_id = hdf5read(edfile,'/PACO_ID');
    paco_n  = hdf5read(edfile,'/PACO_N');
    agb_co  = hdf5read(edfile,'/AGB_CO');
    np_co   = hdf5read(edfile,'/NPLANT');
    ba_co   = hdf5read(edfile,'/BA_CO');
    dbh_co  = hdf5read(edfile,'/DBH');
    pft_co  = hdf5read(edfile,'/PFT');
    dagb_co = hdf5read(edfile,'/DAGB_DT');
    ddbh_co = hdf5read(edfile,'/DDBH_DT');
    bdead_co = hdf5read(edfile,'/BDEAD');
    balive_co = hdf5read(edfile,'/BALIVE');
    dist_type  = hdf5read(edfile,'/DIST_TYPE');
    patch_area = hdf5read(edfile,'/AREA');
    patch_age  = hdf5read(edfile,'/AGE');
    site_area  = hdf5read(edfile,'/AREA_SI');
    
    mmean_mort_rate  = hdf5read(edfile,'/MMEAN_MORT_RATE_CO');
    isi_b = pysi_id(1);
    isi_e = isi_b+pysi_n(1) - 1;
    ba_py  = 0;
    agb_py = 0;
    area_py = 0;
    ip=0;
    
    for isi=isi_b:isi_e
        ipa_b = sipa_id(isi);
        ipa_e = ipa_b + sipa_n(isi) - 1;
        for ipa=ipa_b:ipa_e
            if dist_type(ipa)==3
                ip=ip+1;
                mdata.date(iy,im).pa(ip).afrac = (site_area(isi).*patch_area(ipa))/areaq(it);
                mdata.date(iy,im).pa(ip).age   = patch_age(ipa);
                
                % HERE WE ARE ALSO FILTERING OUT THE <10cm FELLOWS
                ico_b = int16(paco_id(ipa));
                ico_e = ico_b + int16(paco_n(ipa)) - 1;
                
                local_dbs = int16(find(dbh_co(ico_b:ico_e)>dbhmin));
                ico_dbs = ico_b+local_dbs-1;
                
                mdata.date(iy,im).pa(ip).ed_con    = np_co(ico_dbs);
                mdata.date(iy,im).pa(ip).ed_copft  = pft_co(ico_dbs);
                mdata.date(iy,im).pa(ip).ed_codbh  = dbh_co(ico_dbs);
                mdata.date(iy,im).pa(ip).ed_coba   = ba_co(ico_dbs);
                mdata.date(iy,im).pa(ip).ed_coagb  = agb_co(ico_dbs);
                mdata.date(iy,im).pa(ip).ed_cobm   = bdead_co(ico_dbs)+balive_co(ico_dbs);
                mdata.date(iy,im).pa(ip).ed_dagb   = dagb_co(ico_dbs);
                mdata.date(iy,im).pa(ip).ed_ddbh   = ddbh_co(ico_dbs);
                mdata.date(iy,im).pa(ip).ed_mr     = mmean_mort_rate(:,ico_dbs);
                
                min_npaf = double(mdata.date(iy,im).pa(ip).afrac)*double(min(mdata.date(iy,im).pa(ip).ed_con));
                if(min_npaf<mdata.date(iy,im).min_npaf_all)
                    mdata.date(iy,im).min_npaf_all=min_npaf;
                end
                
                for ico=ico_dbs'
                    if (pft_co(ico)==2 && dbh_co(ico)>max_dbh_early)
                        max_dbh_early=dbh_co(ico);
                    end
                    if (pft_co(ico)==4 && dbh_co(ico)>max_dbh_late)
                        max_dbh_late=dbh_co(ico);
                    end
                    
                end
                
            end % if dist_type
        end % ipa
    end %isi
    
    mdata.date(iy,im).npatch = ip;
    
end

display('Finished loading model Q Files');
