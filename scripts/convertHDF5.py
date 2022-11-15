import numpy as np
import h5py
import os
import glob

for filename in glob.glob('data/QP_BSE/ndb.BS*'):

  if("hdf5" in filename or ("PAR_Q" not in filename and "head_Q1" not in filename)):
    continue

  # creating a file
  oldFile = h5py.File(filename, 'r')
  newFile = h5py.File(filename+'.hdf5','w')

  print(filename)

  # first header file contains different info
  if("head_Q1" in filename):

    for ik, key in enumerate(['HEAD_QPT','Bands']): #fil.keys()):

      # open the data
      data = np.array(list(oldFile[key]))
      if(key == 'Bands'):
        data = data.astype(np.int)
        data = np.array([[data[0],data[1]]])
        size = (2)
      #if(key == 'HEAD_QPT'):
      #  data = data.astype(np.float64)
      #  size = data.shape
      newFile.create_dataset(key, data=data);

  elif ("PAR_Q" in filename):

    for ik, key in enumerate(['BSE_RESONANT']):

      # open the data
      data = np.array(list(oldFile[key]))
      if(key == 'BSE_RESONANT'):
        data = 1j*data[...,1] + data[...,0]
      newFile.create_dataset(key, data=data);

  newFile.close()
  oldFile.close()


