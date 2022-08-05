import xarray as xr
import numpy as np

lat_pt = 88.3935
lon_pt = 80.51758

mems = np.arange(1,81)


data = xr.open_dataset('/glade/p/cesmdata/cseg/inputdata/ocn/docn7/SOM/pop_frc.b.e21.BW1850.f09_g17.CMIP6-piControl.001.190514.nc')
lon = data.xc.values
lat = data.yc.values
a = abs(lon-lon_pt) + abs(lat-lat_pt)
i,j = np.unravel_index(a.argmin(),a.shape)
T = data.T.values[:,i,j]
S = data.S.values[:,i,j]
hblt = data.hblt.values[:,i,j]
U = data.U.values[:,i,j]
V = data.V.values[:,i,j]
dhdx = data.dhdx.values[:,i,j]
dhdy = data.dhdy.values[:,i,j]
qdp = data.qdp.values[:,i,j]
data.close()


for m in range(0,mems.shape[0]):
  label = str(m+1).zfill(4)
  print('Working on => {0}'.format(label))
  new_file = open('OCN_FORCING_{0}.txt'.format(label),'w')
  new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(T[0],T[1],T[2],T[3],T[4],T[5],T[6],T[7],T[8],T[9],T[10],T[11]))
  new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(S[0],S[1],S[2],S[3],S[4],S[5],S[6],S[7],S[8],S[9],S[10],S[11]))
  new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(hblt[0],hblt[1],hblt[2],hblt[3],hblt[4],hblt[5],hblt[6],hblt[7],hblt[8],hblt[9],hblt[10],hblt[11]))
  new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(U[0],U[1],U[2],U[3],U[4],U[5],U[6],U[7],U[8],U[9],U[10],U[11]))
  new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(V[0],V[1],V[2],V[3],V[4],V[5],V[6],V[7],V[8],V[9],V[10],V[11]))
  new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(dhdx[0],dhdx[1],dhdx[2],dhdx[3],dhdx[4],dhdx[5],dhdx[6],dhdx[7],dhdx[8],dhdx[9],dhdx[10],dhdx[11]))
  new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(dhdy[0],dhdy[1],dhdy[2],dhdy[3],dhdy[4],dhdy[5],dhdy[6],dhdy[7],dhdy[8],dhdy[9],dhdy[10],dhdy[11]))
  new_file.writelines('{0:10.13f} {1:10.13f} {2:10.13f}\n    {3:10.13f} {4:10.13f} {5:10.13f}\n    {6:10.13f} {7:10.13f} {8:10.13f}\n    {9:10.13f} {10:10.13f} {11:10.13f}\n'.format(qdp[0],qdp[1],qdp[2],qdp[3],qdp[4],qdp[5],qdp[6],qdp[7],qdp[8],qdp[9],qdp[10],qdp[11]))
  new_file.close() 

