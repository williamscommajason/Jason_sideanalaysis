import rice_decode
import scipy.fftpack as fft
import numpy as np
import scipy.io as io

def dct_decode(dct_error,indices,values):

    
    if dct_error[0] == 1:

        dct_error.pop(0)
        for i in range(len(dct_error))[1:]:
            dct_error[i] += dct_error[i-1]
                
    else:
        dct_error.pop(0)
        
    if values[0] == 1:

        values.pop(0)
        for i in range(len(values))[1:]:
            values[i] += values[i-1]

    else:
        values.pop(0)


    if indices[0] == 1:
        print(indices)
        indices.pop(0)
        for i in range(len(indices))[1:]:
            indices[i] += indices[i-1]

    else:
        indices.pop(0)
        

    xerror = np.zeros(len(dct_error))

    for i in range(len(indices)):
        xerror[int(indices[i])] = values[int(i)]

    xx = fft.idct(xerror,norm='ortho')
    
    error = -(dct_error - xx)
    error = [int(round(x)) for x in error]
    
   
    return error


if __name__ == '__main__':

    error = dct_decode()
    oerror = io.loadmat('error.mat')['error'][0]
    oerror = [int(x) for x in oerror]
#    print((np.array(oerror)-np.array(error)).tolist())

