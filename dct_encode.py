import scipy.fftpack as fft
import numpy as np 
from numpy.linalg import norm
import rice_encode 

def get_dct(x):
    
    xx = fft.dct(x,norm='ortho')
    dx = np.zeros((len(x),2))
    dx[:,1] = fft.dct(x,norm='ortho')
    dx[:,0] = list(range(len(dx)))
    dx = dx[abs(dx[:,1]).argsort()]
    dx = np.flip(dx,0)

    i = 0
    
    while norm(dx[:i,1])/norm(dx[:,1]) < .99:

        i += 1
        norm(dx[:i,1])/norm(dx[:i,1])

    needed = i
     
  
    return dx[0:needed+2,0],dx[0:needed+2,1]

def dct_encode(x):

    indices,values = get_dct(x)
    
    indices = [int(x) for x in indices]
    values = [int(x) for x in values]

    xerror = np.zeros(len(x))
    
    for i in range(len(indices)):
        xerror[int(indices[i])] = values[int(i)]

    xx = fft.idct(xerror,norm='ortho')
    error = [int(round(j)) for j in xx-x]
   

    if np.var(error) < np.var(np.diff(error)):
        error.insert(0,0) 
        #rice_encode.compress(error,'dct_error.bin')
    else:
        derror = [int(x) for x in np.diff(error)]
        derror.insert(0,1)
        derror.insert(1,error[0])
        error = derror
        #rice_encode.compress(derror,'dct_error.bin')

    if np.var(indices) < np.var(np.diff(indices)):
        indices.insert(0,0)
        #rice_encode.compress(indices,'indices.bin')
    else:
        dindices = [int(x) for x in np.diff(indices)]
        dindices.insert(0,1)
        dindices.insert(1,indices[0])
        indices = dindices
        #rice_encode.compress(dindices,'indices.bin')

    if np.var(values) < np.var(np.diff(values)):
        values.insert(0,0)
        #rice_encode.compress(values,'values.bin')
    else:
        dvalues = [int(x) for x in np.diff(values)]
        dvalues.insert(0,1)
        dvalues.insert(1,values[0])
        values = dvalues
        #rice_encode.compress(dvalues,'values.bin')


    return error, indices, values

if __name__ == '__main__':

    import scipy.io as io
    error = io.loadmat('error.mat')['error'][0]
    error = [int(x) for x in error]
    dct_error,indices,values = dct_encode(error)
    f = open('error' + '.emd.dct', 'wb')
    #f = rice_encode.compress(emd.error,f)
    f = rice_encode.compress(dct_error,f)
    f = rice_encode.compress(indices, f)
    f = rice_encode.compress(values,f)
    f.close()
    
