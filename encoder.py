from LPC import LPC
from EMD import EMD
import rice_encode
import rice_decode
import argparse
import numpy as np
import math
import dct_encode
import os
np.set_printoptions(threshold=np.nan)

def encoder():

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()
       
    x = np.load(args.filename)
    emd = run_EMD_encoder(x,len(x))
    r_err = run_LPC_encoder(emd.residual,args.filename)
    emd.make_lossless(r_err)
    dct_error,indices,values = dct_encode.dct_encode(emd.error)

    with open(args.filename + '.emd.dct', 'wb') as f:
        f = rice_encode.compress(dct_error,f)
        f = rice_encode.compress(indices, f)
        f = rice_encode.compress(values,f)
    f.close()
    
    with open(args.filename + '.emd', 'wb') as f:
        f = rice_encode.compress(emd.error,f)

    f.close()
    if os.path.getsize('%s' %args.filename + '.emd') > os.path.getsize(args.filename + '.emd.dct'):
       os.remove(args.filename + '.emd')
       filesize = (os.path.getsize(args.filename + '.emd.dct')) + 40
       
    else:
        filesize = os.path.getsize(args.filename + '.emd') + 40
        os.remove(args.filename + '.emd.dct')        
    return filesize    

def run_EMD_encoder(x,samples):

    emd = EMD(samples)
    emd.emd(x,None)
    return emd

def run_LPC_encoder(residual,filename):

    lpc = LPC(2,len(residual),len(residual))
    lpc.lpc_fit(residual)
    lpc.get_fits(lpc.err)
    lpc.get_amp(lpc.err,lpc.h)
    lpc.pack_residual(filename)
    npts,frame_width,amp,gains,fits = lpc.unpack_residual(filename)
    lpc.recon_err(npts,frame_width,amp,fits)
    r_err = lpc.lpc_synth(lpc.aaa,gains,LPC.r_err,npts,frame_width)
    

    return r_err
    
       

 
if __name__ == "__main__":
    print(encoder())
       
