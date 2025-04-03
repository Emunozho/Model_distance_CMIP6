#!/usr/bin/env python
# coding: utf-8

# In[1]:
# All outputs in kgC m-2 y-1

import numpy as np
np.float_ = np.float64
import glob
import os
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import nc_time_axis
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import xesmf as xe
import cartopy.crs as ccrs
import cftime
import gc
from datetime import datetime
gc.collect()
maskls=xr.open_mfdataset('/home/estefania/cmip6stop/land_sea_masks/landsea.nc')
maskl=maskls.LSMASK.where(maskls.LSMASK!=0)==1; maskl=maskl.where(maskl==1)

def readCMIP6data(path_to_files):
    "Takes NetCDF files from a directory and combines them together by coordinates"
    dataset=[]
    for file in os.listdir(path_to_files):
            if file.endswith(".nc"):
                p=os.path.join(path_to_files, file)
                dataset.append(p)
    out=xr.open_mfdataset(dataset, combine='by_coords')
    return out;

def to_360day_monthly(da):
    '''Takes a DataArray. Change the 
    calendar to 360_day and precision to monthly.'''
    val = da.copy()
    time1 = da.time.copy()
    for itime in range(val.sizes['time']):
        bb = val.time.values[itime].timetuple()
        time1.values[itime] = cftime.Datetime360Day(bb[0],bb[1],1)

    # We rename the time dimension and coordinate to time360 to make it clear it isn't 
    # the original time coordinate.
    val = val.rename({'time':'time'})
    time1 = time1.rename({'time':'time'})
    val = val.assign_coords({'time':time1})
    return val
    
def add_time_dim(xda):
    xda = xda.expand_dims(time = [datetime.now()])
    return xda

def get_nepCMIP6(d1,d2):
    institution=['BCC','CCCma','CNRM-CERFACS','CSIRO','MOHC','MPI-M','NCAR','NCC','NOAA-GFDL']
    model=['BCC-CSM2-MR','CanESM5','CNRM-ESM2-1','ACCESS-ESM1-5','UKESM1-0-LL','MPI-ESM1-2-LR','CESM2','NorESM2-LM','GFDL-ESM4']
    model2=['BCC','CCCma','CNRM','CSIRO','MOHC','MPI','NCAR','NCC','NOAA']

    nf=len(model)
# Land-sea masks
    pth_masks=glob.glob('/home/estefania/cmip6stop/land_sea_masks/**') # read sea-land masks.
    namemask=[]
    for j in range(1,len(pth_masks)):namemask.append(pth_masks[j].split('_')[4])
    
    for i in range(nf):
     #   print('Model:'+model[i])
# CESM2 does not have gpp, rh, ra and mrso.
        if len(glob.glob('/home/data/CMIP6/data/CMIP/'+institution[i]+'/'+model[i]+'/esm-hist/**/Lmon/gpp/**/**'))!=0:
        #if len(glob.glob('/home/data/CMIP6/data/CMIP/'+institution[i]+'/'+model[i]+'/esm-hist/**/Lmon/gpp/**/**/**')>1)
            pth_gpp=glob.glob('/home/data/CMIP6/data/CMIP/'+institution[i]+'/'+model[i]+'/esm-hist/**/Lmon/gpp/**/**')
            pth_gpp=pth_gpp[len(pth_gpp)-1] if len(pth_gpp)>1 else pth_gpp[0]# in some models there are more than 1 version.
            pth_ra=glob.glob('/home/data/CMIP6/data/CMIP/'+institution[i]+'/'+model[i]+'/esm-hist/**/Lmon/ra/**/**')
            pth_ra=pth_ra[len(pth_ra)-1] if len(pth_ra)>1 else pth_ra[0]
            pth_rh=glob.glob('/home/data/CMIP6/data/CMIP/'+institution[i]+'/'+model[i]+'/esm-hist/**/Lmon/rh/**/**')
            pth_rh=pth_rh[len(pth_rh)-1] if len(pth_rh)>1 else pth_rh[0]
            nep_ind=1
        else:nep_ind=0
        if nep_ind==1:
            gpp=readCMIP6data(pth_gpp) # gross primary productivity of carbon [kgC m-2 s-1]. 
            gpp=gpp.where(gpp.gpp!=1e+20)#; gpp=gpp.where(gpp.gpp>=0)
            ra=readCMIP6data(pth_ra) # autotrophic respiration carbon flux [kgC m-2 s-1].
            ra=ra.where(ra.ra!=1e+20)#; ra=ra.where(ra.ra>=0)
            rh=readCMIP6data(pth_rh) # heterotrophic_respiration_carbon_flux [kgC m-2 s-1]. 
            rh=rh.where(rh.rh!=1e+20)#; rh=rh.where(rh.rh>=0)
# Select only the period in common with the product based on observations.
            gpp=gpp.sel(time=slice(d1,d2)) # common period FLUXCOM-X.
            ra=ra.sel(time=slice(d1,d2)) # common period FLUXCOM-X.
            rh=rh.sel(time=slice(d1,d2)) # common period FLUXCOM-X.
# mask the land. 
            if model[i] in namemask:
                pth_mask=pth_masks[namemask.index(model[i])+1]
                mask=xr.open_mfdataset(pth_mask,combine='by_coords')
                if i==7: mask=mask.reindex_like(gpp, method='nearest',tolerance=0.01) # because the mask does not have exactly the same coordinates.
                mask=mask.sftlf.where(mask.sftlf>5)*0+1
                if nep_ind==1:
                    gpp=gpp*mask
                    ra=ra*mask
                    rh=rh*mask
### Regridding to 1$^\circ$ x 1$^\circ$
            ds_target=maskls.LSMASK
            regridder=xe.Regridder(gpp.gpp,ds_target,"bilinear",periodic=True)  
            gpp_reg=regridder(gpp.gpp) 
            ra_reg=regridder(ra.ra)  
            rh_reg=regridder(rh.rh) 
# Change of units
            gpp_orig=gpp; gpp=gpp_reg*(60*60*24*364.3) # Change units to kgC m-2 yr-1.
            ra_orig=ra; ra=ra_reg*(60*60*24*364.3) # Change units to kgC m-2 yr-1.
            rh_orig=rh; rh=rh_reg*(60*60*24*364.3) # Change units to kgC m-2 yr-1.
# Same calendar
            gpp=gpp.convert_calendar("noleap",align_on='year')
            rh=rh.convert_calendar("noleap",align_on='year')
            ra=ra.convert_calendar("noleap",align_on='year')        
# Align to the same day of the month.
            gpp=to_360day_monthly(gpp)
            rh=to_360day_monthly(rh)
            ra=to_360day_monthly(ra)
# NEP
            nep=gpp-(ra+rh) # kgC m-2 yr-1
            if i==0:
                nep_cmip6=xr.Dataset({model2[0]:nep})
                gpp_cmip6=xr.Dataset({model2[0]:gpp})
                ra_cmip6=xr.Dataset({model2[0]:ra})
                rh_cmip6=xr.Dataset({model2[0]:rh})
            else:
                nep_cmip6=xr.merge([nep_cmip6,nep.to_dataset(name=model2[i])])
                gpp_cmip6=xr.merge([gpp_cmip6,gpp.to_dataset(name=model2[i])])
                ra_cmip6=xr.merge([ra_cmip6,ra.to_dataset(name=model2[i])])
                rh_cmip6=xr.merge([rh_cmip6,rh.to_dataset(name=model2[i])])
    return(nep_cmip6,gpp_cmip6,ra_cmip6,rh_cmip6)

def get_nbpCMIP6(d1,d2):
    institution=['BCC','CCCma','CNRM-CERFACS','CSIRO','MOHC','MPI-M','NCAR','NCC','NOAA-GFDL']
    model=['BCC-CSM2-MR','CanESM5','CNRM-ESM2-1','ACCESS-ESM1-5','UKESM1-0-LL','MPI-ESM1-2-LR','CESM2','NorESM2-LM','GFDL-ESM4']
    model2=['BCC','CCCma','CNRM','CSIRO','MOHC','MPI','NCAR','NCC','NOAA']

    nf=len(model)
# Land-sea masks
    pth_masks=glob.glob('/home/estefania/cmip6stop/land_sea_masks/**') # read sea-land masks.
    namemask=[]
    for j in range(1,len(pth_masks)):namemask.append(pth_masks[j].split('_')[4])
    
    for i in range(nf):
        if len(glob.glob('/home/data/CMIP6/data/CMIP/'+institution[i]+'/'+model[i]+'/esm-hist/**/Lmon/nbp/**/**'))!=0:
            pth_nbp=glob.glob('/home/data/CMIP6/data/CMIP/'+institution[i]+'/'+model[i]+'/esm-hist/**/Lmon/nbp/**/**')
            pth_nbp=pth_nbp[len(pth_nbp)-1] if len(pth_nbp)>1 else pth_nbp[0]
            nbp_ind=1
        else:nbp_ind=0
        
        if nbp_ind==1:
            nbp=readCMIP6data(pth_nbp)
            nbp=nbp.where(nbp.nbp!=1e+20)
            nbp=nbp.sel(time=slice(d1,d2)) # common period CarboScope.
# mask the land. 
            if model[i] in namemask:
                pth_mask=pth_masks[namemask.index(model[i])+1]
                mask=xr.open_mfdataset(pth_mask,combine='by_coords')
                if i==7: mask=mask.reindex_like(nbp,method='nearest',tolerance=0.01) # because the mask does not have exactly the same coordinates.
                mask=mask.sftlf.where(mask.sftlf>5)*0+1
                nbp=nbp*mask
            ds_target=maskls.LSMASK
            regridder=xe.Regridder(nbp.nbp,ds_target,"bilinear",periodic=True) 
            nbp_reg=regridder(nbp.nbp) 
            nbp_orig=nbp; nbp=nbp_reg*(60*60*24*364.3) # Change units to kgC m-2 yr-1.
            nbp=nbp.convert_calendar("noleap",align_on='year')
            nbp=to_360day_monthly(nbp)
            if i==1:nbp_cmip6=xr.Dataset({model2[1]:nbp})
            else:nbp_cmip6=xr.merge([nbp_cmip6,nbp.to_dataset(name=model2[i])])
    return(nbp_cmip6)

def get_neeFLUXCOM(d1,d2):
    pth=glob.glob('/home/data/Fluxcom/NEE/NEE_????_005_monthly.nc')
    pth.sort()
    nee=xr.open_mfdataset(pth,combine='by_coords') # [gC m-2 d-1].
    nee=nee.sel(time=slice(d1,d2)) # Select only the period in common with the product based on modeling.
    ds_target=maskls.LSMASK
    regridder=xe.Regridder(nee.NEE,ds_target,"bilinear",periodic=True)
    nee_reg=regridder(nee.NEE) 
    nee=nee_reg*364.3/1e3 # Change units to kgC m-2 yr-1.
    nee=nee.convert_calendar("noleap",align_on='year') # same calendar
    nee_fc=to_360day_monthly(nee) # Align to the same day of the month
    return(nee_fc)

def get_nbpCarboscope(d1,d2):
    df=xr.open_dataset('/home/data/Inversion_models/CarboScope/nbetEXToc_v2024E.flux.nc') #  [PgC/yr]
    nbp=df.co2flux_land/df.dxyp #  [PgC/yr/m2]
    nbp=nbp.rename(mtime="time")
    nbp=nbp.sel(time=slice(d1,d2)) # Select only the period in common with the product based on modeling.
    ds_target=maskls.LSMASK
    regridder=xe.Regridder(nbp,ds_target,"bilinear",periodic=True)
    nbp_reg=regridder(nbp) 
    nbp=nbp_reg*1e+12 # Change units to kg C m-2 yr-1.
    nbp=nbp.resample(time="ME",label="right").mean(skipna=False,keep_attrs=True).compute()
    nbp=nbp*maskl
    nbp=nbp.convert_calendar("noleap",align_on='year')
    nbp_CaS=to_360day_monthly(nbp) # Align to the same day of the month
    return(nbp_CaS)

def get_neeCarboscope(d1,d2):
    df=xr.open_zarr('/home/data/Inversion_models/nee_050_monthly',consolidated=True) # This only works if xarray[complete] is installed.
    nee=df.CarboScope # [gC m-2 d-1] Here it is assumed that the value given by CarboScope minus the contributions from fire equals NEE.
    ds_target=maskls.LSMASK
    regridder=xe.Regridder(nee,ds_target,"bilinear",periodic=True)
    nee_reg=regridder(nee) 
    nee=nee_reg*364.3/1e+3 # Change units to kg C m-2 yr-1.
    nee=nee.sel(time=slice(d1,d2)) # Select only the period in common with the product based on modeling.
    nee=nee.convert_calendar("noleap",align_on='year') # same calendar
    nee_CaS=to_360day_monthly(nee) 
    return(nee_CaS)

def get_nbpCAMS(d1,d2):
    pth=glob.glob('/home/data/Inversion_models/CAMS/CAMSv24/*.nc')
    pth.sort()
    dates=pd.date_range(datetime.strptime("1979-01-01","%Y-%m-%d"),datetime.strptime("2020-01-01","%Y-%m-%d"),freq='ME')
    nbp=xr.open_mfdataset(pth,preprocess=add_time_dim)
    nbp['time']=dates
    nbp=nbp.flux_apos_bio #[kgC m-2 month-1]
    nbp=nbp.sel(time=slice(d1,d2)) # Select only the period in common with the product based on modeling.
    ds_target=maskls.LSMASK
    regridder=xe.Regridder(nbp,ds_target,"bilinear",periodic=True)
    nbp_reg=regridder(nbp) 
    nbp_orig=nbp 
    nbp=nbp_reg*12 # Change units to kgC m-2 yr-1 
    nbp=nbp*maskl
    nbp=nbp.convert_calendar("noleap",align_on='year')
    nbp_CAMS=to_360day_monthly(nbp) # Align to the same day of the month
    return(nbp_CAMS)

def get_nepTrendy(d1,d2): 
    model=['CABLE-POP','CARDAMOM','CLASSIC','CLM5.0','DLEM','ED','ELM','IBIS','iMAPLE','ISAM','ISBACTRIP','JSBACH','JULES','LPJml','LPJwsl','LPX','OCN','ORCHIDEE','SDGVM','VISIT','VISIT-UT']
    nf=len(model)
    for i in range(nf): #nf
        print(model[i])
        gc.collect()
        pth_gpp=glob.glob('/home/data/Trendy/gcb2024/'+model[i]+'/*gpp.nc')[0]
        pth_ra=glob.glob('/home/data/Trendy/gcb2024/'+model[i]+'/*ra.nc')[0]
        pth_rh=glob.glob('/home/data/Trendy/gcb2024/'+model[i]+'/*rh.nc')[0]
        if model[i] in ['LPJml','VISIT-UT','CABLE-POP','CARDAMOM','DLEM','ISAM','VISIT','iMAPLE']:
            gpp=xr.open_dataset(pth_gpp,decode_times=False)
            ra=xr.open_dataset(pth_ra,decode_times=False)
            rh=xr.open_dataset(pth_rh,decode_times=False)
            if model[i] in ['ISAM','VISIT-UT']:gpp.time.attrs['units'] = 'months since 1700-1-1 12:00:00'
            if model[i]=='LPJml': units,reference_date=gpp.time.attrs['units'].split('since')
            elif model[i]=='CARDAMOM':units,reference_date=gpp.Time.attrs['long_name'].split('since')
            elif model[i]=='DLEM':units,reference_date=gpp.time.attrs['unit'].split('since')
            else:units,reference_date=gpp.time.attrs['units'].split('since')
            gpp['time']=pd.date_range(start=reference_date,periods=gpp.sizes['time'],freq='MS')
            ra['time']=pd.date_range(start=reference_date,periods=ra.sizes['time'],freq='MS')
            rh['time']=pd.date_range(start=reference_date,periods=rh.sizes['time'],freq='MS')       
        else: 
            gpp=xr.open_dataset(pth_gpp)
            ra=xr.open_dataset(pth_ra)
            rh=xr.open_dataset(pth_rh)
        if model[i] in ['CLM5.0','DLEM','ELM','ISAM','ISBACTRIP','JULES','VISIT','VISIT-UT']:
            if model[i] in ['ELM']:
                gpp=gpp.set_index(lat="latitude",lon='longitude')
                ra=ra.set_index(lat="latitude",lon='longitude')
                rh=rh.set_index(lat="latitude",lon='longitude')
            else:
                gpp=gpp.rename(lat="latitude",lon='longitude')
                ra=ra.rename(lat="latitude",lon='longitude')
                rh=rh.rename(lat="latitude",lon='longitude') 
# Select only the period in common.
        gpp=gpp.sel(time=slice(d1,d2))
        ra=ra.sel(time=slice(d1,d2))
        rh=rh.sel(time=slice(d1,d2))
### Regridding to 1$^\circ$ x 1$^\circ$
        ds_target=maskls.LSMASK
        regridder=xe.Regridder(gpp.gpp,ds_target,"bilinear",periodic=True)  
        gpp_reg=regridder(gpp.gpp) 
        ra_reg=regridder(ra.ra)  
        rh_reg=regridder(rh.rh) 
# Change of units
        gpp_orig=gpp; gpp=gpp_reg*(60*60*24*364.3)*maskl # Change units to kgC m-2 yr-1. 
        ra_orig=ra; ra=ra_reg*(60*60*24*364.3)*maskl # Change units to kgC m-2 yr-1.
        rh_orig=rh; rh=rh_reg*(60*60*24*364.3)*maskl # Change units to kgC m-2 yr-1.
# # Same calendar
        gpp=gpp.convert_calendar("noleap",align_on='year')
        rh=rh.convert_calendar("noleap",align_on='year')
        ra=ra.convert_calendar("noleap",align_on='year')        
# Align to the same day of the month.
        gpp=to_360day_monthly(gpp)
        rh=to_360day_monthly(rh)
        ra=to_360day_monthly(ra)
# # NEP
        nep=gpp-(ra+rh) # kgC m-2 yr-1
# Merge all models.
        if i==0:
            nep_trendy=xr.Dataset({model[0]:nep})
            gpp_trendy=xr.Dataset({model[0]:gpp})
            ra_trendy=xr.Dataset({model[0]:ra})
            rh_trendy=xr.Dataset({model[0]:rh})
        else:
            nep_trendy=xr.merge([nep_trendy,nep.to_dataset(name=model[i])])
            gpp_trendy=xr.merge([gpp_trendy,gpp.to_dataset(name=model[i])])
            ra_trendy=xr.merge([ra_trendy,ra.to_dataset(name=model[i])])
            rh_trendy=xr.merge([rh_trendy,rh.to_dataset(name=model[i])])
    return(nep_trendy,gpp_trendy,ra_trendy,rh_trendy)



    
#     # ['CABLEPOP','CLASSIC','CLM5.0','DLEM','ED','ELM','IBIS','ISAM','ISBACTRIP','JSBACH','JULES','LPJml','LPJwsl','OCN','ORCHIDEE','SDGVM','VISIT','YIBS','lpxqs']
                           
# # 'LPJ-GUESS' does not have rh. CARDAMOM only since 2003.
#     nf=len(model)
#     for i in range(nf): #nf
#         #print(model[i])
#         pth_gpp=glob.glob('/home/data/Trendy/Trendy_v13/'+model[i]+'/S3/*gpp*.nc')[0]
#         pth_ra=glob.glob('/home/data/Trendy/Trendy_v13/'+model[i]+'/S3/*ra.nc')[0]
#         pth_rh=glob.glob('/home/data/Trendy/Trendy_v13/'+model[i]+'/S3/*rh*.nc')[0]
#         if model[i] in ['CABLEPOP','CARDAMOM','DLEM','IBIS','ISAM','ISBACTRIP','LPJml','VISIT','YIBS']:
#             gpp=xr.open_dataset(pth_gpp,decode_times=False)
#             ra=xr.open_dataset(pth_ra,decode_times=False)
#             rh=xr.open_dataset(pth_rh,decode_times=False)
#             if model[i]=='ISBACTRIP':
#                 gpp=gpp.rename(lat_FULL="latitude",lon_FULL='longitude',time_counter='time')
#                 ra=ra.rename(lat_FULL="latitude",lon_FULL='longitude',time_counter='time')
#                 rh=rh.rename(lat_FULL="latitude",lon_FULL='longitude',time_counter='time')
#             elif model[i]=='CARDAMOM':
#                 units,reference_date=gpp.Time.attrs['long_name'].split('since')
# #            gpp=gpp.rename(Time='time');ra=ra.rename(Time='time');rh=rh.rename(Time='time')
#             elif model[i]=='DLEM':units,reference_date=gpp.time.attrs['unit'].split('since')
#             else:units,reference_date=gpp.time.attrs['units'].split('since')
#             gpp['time']=pd.date_range(start=reference_date,periods=gpp.sizes['time'],freq='MS')
#             ra['time']=pd.date_range(start=reference_date,periods=ra.sizes['time'],freq='MS')
#             rh['time']=pd.date_range(start=reference_date,periods=rh.sizes['time'],freq='MS')     
#         else:
#             gpp=xr.open_dataset(pth_gpp)
#             ra=xr.open_dataset(pth_ra)
#             rh=xr.open_dataset(pth_rh)
#         if model[i] not in ['CLM5.0','DLEM','ISAM','JULES','LPJwsl','VISIT']:
#             if model[i] in ['ELM']:
#                 gpp=gpp.set_index(lat="latitude",lon='longitude')
#                 ra=ra.set_index(lat="latitude",lon='longitude')
#                 rh=rh.set_index(lat="latitude",lon='longitude')
#             else:
#                 gpp=gpp.rename(latitude="lat",longitude='lon')
#                 ra=ra.rename(latitude="lat",longitude='lon')
#                 rh=rh.rename(latitude="lat",longitude='lon') 
# # Select only the period in common.
#         gpp=gpp.sel(time=slice(d1,d2))
#         ra=ra.sel(time=slice(d1,d2))
#         rh=rh.sel(time=slice(d1,d2))
# ### Regridding to 1$^\circ$ x 1$^\circ$
#         ds_target=maskls.LSMASK
#         regridder=xe.Regridder(gpp.gpp,ds_target,"bilinear",periodic=True)  
#         gpp_reg=regridder(gpp.gpp) 
#         ra_reg=regridder(ra.ra)  
#         rh_reg=regridder(rh.rh) 
# # Change of units
#         gpp_orig=gpp; gpp=gpp_reg*(60*60*24*364.3)*maskl # Change units to kgC m-2 yr-1. 
#         ra_orig=ra; ra=ra_reg*(60*60*24*364.3)*maskl # Change units to kgC m-2 yr-1.
#         rh_orig=rh; rh=rh_reg*(60*60*24*364.3)*maskl # Change units to kgC m-2 yr-1.
# # # Same calendar
#         gpp=gpp.convert_calendar("noleap",align_on='year')
#         rh=rh.convert_calendar("noleap",align_on='year')
#         ra=ra.convert_calendar("noleap",align_on='year')        
# # Align to the same day of the month.
#         gpp=to_360day_monthly(gpp)
#         rh=to_360day_monthly(rh)
#         ra=to_360day_monthly(ra)
# # # NEP
#         nep=gpp-(ra+rh) # kgC m-2 yr-1
# # Merge all models.
#         if i==0:
#             nep_trendy=xr.Dataset({model[0]:nep})
#             gpp_trendy=xr.Dataset({model[0]:gpp})
#             ra_trendy=xr.Dataset({model[0]:ra})
#             rh_trendy=xr.Dataset({model[0]:rh})
#         else:
#             nep_trendy=xr.merge([nep_trendy,nep.to_dataset(name=model[i])])
#             gpp_trendy=xr.merge([gpp_trendy,gpp.to_dataset(name=model[i])])
#             ra_trendy=xr.merge([ra_trendy,ra.to_dataset(name=model[i])])
#             rh_trendy=xr.merge([rh_trendy,rh.to_dataset(name=model[i])])
#     return(nep_trendy,gpp_trendy,ra_trendy,rh_trendy)
    
def get_nbpTrendy(d1,d2): 
    model=['CABLE-POP','CARDAMOM','CLASSIC','CLM5.0','ED','ELM','IBIS','iMAPLE','ISAM','ISBACTRIP','JSBACH','JULES','LPJwsl','LPX','OCN','ORCHIDEE','SDGVM','VISIT','VISIT-UT']
    
#    model=['CABLEPOP','CLASSIC','CLM5.0','ED','ELM','IBIS','ISAM','ISBACTRIP','JSBACH','JULES','LPJwsl','OCN','ORCHIDEE','VISIT','YIBS','lpxqs']
# LPJml, DLEM, and LPJ-GUESS only have annual data. CARDADOM only since 2003.
    nf=len(model)
    for i in range(nf): #nf
        print(model[i])
        pth_nbp=glob.glob('/home/data/Trendy/gcb2024/'+model[i]+'/*nbp.nc')[0]

        if model[i] in ['LPJml','VISIT-UT','CABLE-POP','CARDAMOM','DLEM','ISAM','VISIT','iMAPLE']:
            nbp=xr.open_dataset(pth_nbp,decode_times=False)
            if model[i] in ['ISAM','VISIT-UT']:nbp.time.attrs['units']='months since 1700-1-1 12:00:00'
            if model[i]=='LPJml': units,reference_date=nbp.time.attrs['units'].split('since')
            elif model[i]=='CARDAMOM':units,reference_date=nbp.Time.attrs['long_name'].split('since')
            elif model[i]=='DLEM':units,reference_date=nbp.time.attrs['unit'].split('since')
            else:units,reference_date=nbp.time.attrs['units'].split('since')
            nbp['time']=pd.date_range(start=reference_date,periods=nbp.sizes['time'],freq='MS')  
        else: nbp=xr.open_dataset(pth_nbp)
        if model[i] in ['CLM5.0','DLEM','ELM','ISAM','ISBACTRIP','JULES','VISIT','VISIT-UT']:
            nbp=nbp.rename(lat="latitude",lon='longitude')
            if model[i] in ['ELM','IBIS']:nbp=nbp.set_index(latitude="latitude",longitude='longitude')
# Select only the period in common.
        nbp=nbp.sel(time=slice(d1,d2))
### Regridding to 1$^\circ$ x 1$^\circ$
        ds_target=maskls.LSMASK
        regridder=xe.Regridder(nbp.nbp,ds_target,"bilinear",periodic=True)  
        nbp_reg=regridder(nbp.nbp) 
# Change of units
        nbp_orig=nbp; nbp=nbp_reg*(60*60*24*364.3)*maskl # Change units to kgC m-2 yr-1. 
# # Same calendar
        nbp=nbp.convert_calendar("noleap",align_on='year')      
# Align to the same day of the month.
        nbp=to_360day_monthly(nbp)
# Merge all models
        if i==0:nbp_trendy=xr.Dataset({model[0]:nbp})
        else:nbp_trendy=xr.merge([nbp_trendy,nbp.to_dataset(name=model[i])])
    return(nbp_trendy)

        
        
        
        
        
        
        
        
#         if model[i] in ['CABLEPOP','CARDAMOM','DLEM','IBIS','ISAM','ISBACTRIP','VISIT','YIBS']:
#             nbp=xr.open_dataset(pth_nbp,decode_times=False)
#             if model[i]=='ISBACTRIP':nbp=nbp.rename(lat_FULL="latitude",lon_FULL='longitude',time_counter='time')
#             elif model[i]=='CARDAMOM':units,reference_date=nbp.Time.attrs['long_name'].split('since')
#             elif model[i]=='DLEM':units,reference_date=nbp.time.attrs['unit'].split('since')
#             else:units,reference_date=nbp.time.attrs['units'].split('since')
#             nbp['time']=pd.date_range(start=reference_date,periods=nbp.sizes['time'],freq='MS') 
#         else:nbp=xr.open_dataset(pth_nbp)
#         if model[i] not in ['CLM5.0','DLEM','ISAM','JULES','LPJwsl','VISIT']:
#             if model[i] in ['ELM']:nbp=nbp.set_index(lat="latitude",lon='longitude')
#             else: nbp=nbp.rename(latitude="lat",longitude='lon')       
#         if model[i]=='ISAM':nbp=nbp.rename(name_dict={'npp-rh+ld':'nbp'})
# # Select only the period in common.
#         nbp=nbp.sel(time=slice(d1,d2))
# ### Regridding to 1$^\circ$ x 1$^\circ$
#         ds_target=maskls.LSMASK
#         regridder=xe.Regridder(nbp.nbp,ds_target,"bilinear",periodic=True)  
#         nbp_reg=regridder(nbp.nbp) 
# # Change of units
#         nbp_orig=nbp; nbp=nbp_reg*(60*60*24*364.3)*maskl # Change units to kgC m-2 yr-1. 
# # # Same calendar
#         nbp=nbp.convert_calendar("noleap",align_on='year')      
# # Align to the same day of the month.
#         nbp=to_360day_monthly(nbp)
# # Merge all models
#         if i==0:nbp_trendy=xr.Dataset({model[0]:nbp})
#         else:nbp_trendy=xr.merge([nbp_trendy,nbp.to_dataset(name=model[i])])
#    return(nbp_trendy)