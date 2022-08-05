import xarray as xr
import netCDF4
import numpy as np
import glob
import pandas as pd

lat_pt = 88.3935
lon_pt = 80.51758

def convert_times(input_times):
  new_times = np.zeros((input_times.shape[0]),dtype=np.int64)
  for t in range(0,input_times.shape[0]):
    new_times[t] = int(pd.to_datetime(str(input_times[t])).strftime('%Y%m%d%H%M'))
  return new_times

new_times = pd.date_range(start='1/1/2011',end='1/1/2020',freq='1H')[:-1].to_pydatetime()
interp_times = []
for d in new_times:
  if int(d.strftime('%m%d')) == 229:
    continue
  else:
    interp_times.append(int(d.strftime('%Y%m%d%H%M')))
interp_times = np.array(interp_times)



mems = np.arange(1,81)

for m in range(0,mems.shape[0]):
  label = str(mems[m]).zfill(4)
  print('Working on member {0}'.format(label))
  files = xr.open_mfdataset('/glade/collections/rda/data/ds345.0/cpl_unzipped/{0}/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_{0}.ha2x1hi.20*.nc'.format(label))
  if m == 0:
    tlat = files.doma_lat.values[0,:,:]
    tlon = files.doma_lon.values[0,:,:]
    tlon[np.isnan(tlon)] = 0.0
    a = abs(tlat-lat_pt)+abs(tlon-lon_pt)
    i,j = np.unravel_index(a.argmin(),a.shape)
  hi1_times = convert_times(files.time.values)
  hold_sw = files.a2x1hi_Faxa_swndr.values[:,i,j] + files.a2x1hi_Faxa_swvdr.values[:,i,j]+files.a2x1hi_Faxa_swndf.values[:,i,j]+files.a2x1hi_Faxa_swvdf.values[:,i,j]
  files.close()
  files = xr.open_mfdataset('/glade/collections/rda/data/ds345.0/cpl_unzipped/{0}/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_{0}.ha2x1h.20*.nc'.format(label))
  h1_times = convert_times(files.time.values)
  hold_u = files.a2x1h_Sa_u.values[:,i,j] 
  hold_v = files.a2x1h_Sa_v.values[:,i,j]
  files.close()
  files = xr.open_mfdataset('/glade/collections/rda/data/ds345.0/cpl_unzipped/{0}/f.e21.FHIST_BGC.f09_025.CAM6assim.011.cpl_{0}.ha2x3h.20*.nc'.format(label))
  h3_times = convert_times(files.time.values) 
  hold_z = files.a2x3h_Sa_z.values[:,i,j]
  hold_temp = files.a2x3h_Sa_tbot.values[:,i,j]
  hold_dlw = files.a2x3h_Faxa_lwdn.values[:,i,j]
  hold_sphm = files.a2x3h_Sa_shum.values[:,i,j]
  hold_rain = files.a2x3h_Faxa_rainc.values[:,i,j] + files.a2x3h_Faxa_rainl.values[:,i,j]
  hold_snow = files.a2x3h_Faxa_snowc.values[:,i,j] + files.a2x3h_Faxa_snowl.values[:,i,j]
  files.close()


  interp_sw = np.interp(interp_times,hi1_times,hold_sw)
  interp_u = np.interp(interp_times,h1_times,hold_u)
  interp_v = np.interp(interp_times,h1_times,hold_v)
  interp_z = np.interp(interp_times,h3_times,hold_z)
  interp_temp = np.interp(interp_times,h3_times,hold_temp)
  interp_dlw = np.interp(interp_times,h3_times,hold_dlw)
  interp_sphm = np.interp(interp_times,h3_times,hold_sphm)
  interp_rain = np.interp(interp_times,h3_times,hold_rain)
  interp_snow = np.interp(interp_times,h3_times,hold_snow)

  new_file = open('ATM_FORCING_{0}.txt'.format(label),'w')
  new_file.writelines('#Z    DSWSFC      DLWSFC    WNDU10    WNDV10    TEMP2M    SPECHUM    PERCIPRAIN    PRECIPSNOW\n')
  new_file.writelines('# m   W/m**2      w/m**2    m/s       m/s       K         Kg/Kg      kg/m**2/s     kg/m**2/s\n')
  for l in range(0,interp_sw.shape[0]):
    new_file.writelines('{0:10.5f} {1:10.5f} {2:10.5f} {3:10.5f} {4:10.5f} {5:10.5f} {6:10.8f} {7:10.8f} {8:10.8f}\n'.format(interp_z[l],interp_sw[l],interp_dlw[l],interp_u[l],interp_v[l],interp_temp[l],interp_sphm[l],interp_rain[l],interp_snow[l]))
  new_file.close()



