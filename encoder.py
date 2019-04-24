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
    f = open(args.filename + '.emd.dct', 'w+b')
    r_err = run_LPC_encoder(emd.residual,f)
    emd.make_lossless(r_err)
    dct_output = dct_encode.dct_encode(emd.error)

    for i in dct_output:
        f = rice_encode.compress(i,f)
        
    f.close()
    
    with open(args.filename + '.emd', 'w+b') as fd:
        fd = rice_encode.compress(emd.error,fd)

    fd.close()

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

def run_LPC_encoder(residual,f):

    lpc = LPC(2,len(residual),len(residual))
    lpc.lpc_fit(residual)
    lpc.get_fits(lpc.err)
    lpc.get_amp(lpc.err,lpc.h)
    lpc.pack_residual(f)
    npts,frame_width,amp,gains,fits,f = lpc.unpack_residual(f)
    lpc.recon_err(npts,frame_width,amp,fits)
    r_err = lpc.lpc_synth(lpc.aaa,gains,LPC.r_err,npts,frame_width)
    

    return r_err
    
       

 
if __name__ == "__main__":
    print(encoder())
       
