import numpy as np

mems = np.arange(1,81)
truth_mem = 40
ind = np.where(mems==truth_mem)[0][0]
mems = np.delete(mems,ind)

new_file = open('restarts_out.txt','w')
for m in range(0,mems.shape[0]):
  line = str(mems[m]).zfill(4)
  new_file.writelines('mem{0}/restart_state.nc\n'.format(line))
new_file.close()


