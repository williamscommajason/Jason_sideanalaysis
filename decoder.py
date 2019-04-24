from LPC import LPC
import os
import rice_decode
import argparse
import dct_decode

def decoder():

    parser = argparse.ArgumentParser()
    parser.add_argument('filename')
    args = parser.parse_args()

    npts,frame_width,amp,gains,fits = LPC.unpack_residual(args.filename)
    lpc = LPC(2,frame_width,npts)
    LPC.recon_err(npts,frame_width,amp,fits)
    r_err = lpc.lpc_synth(lpc.aaa,gains,LPC.r_err,npts,frame_width)
    if os.path.isfile(args.filename + '.emd.dct'):
        lists = rice_decode.decompress(args.filename + '.emd.dct')
        error = dct_decode.dct_decode(lists[0],lists[1],lists[2])
    
    else:
        error = rice_decode.decompress(args.filename + '.emd')[0]
    
     
    r_sig = error + r_err
   
    return [int(round(x)) for x in r_sig]


if __name__ == "__main__":

   signal = decoder()
   import numpy as np
   np.set_printoptions(threshold=np.nan)
   ts = np.load('timestream1007.npy')
   print(signal-ts)

        
