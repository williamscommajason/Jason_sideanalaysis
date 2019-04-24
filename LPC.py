import logging
import numpy as np
import scipy
import scipy.signal as sig
from scipy.io import loadmat
import math
import struct

class LPC:
    """
    LPC:

    **Linear Predictive Coding**

    Method of predicting a signal based on past samples

    """

    logger = logging.getLogger(__name__)

    def __init__(self, order, frame_width, npts):


        self.aa = [1, -2, 1]
        self.gains = None
        self.h = frame_width
        self.npts = npts
        self.d = None
        self.order = order
        self.amp = None
        self.fits = None
        self.err = None
        self.r_err = None
        self.a = None
        self.nhops = math.floor(self.npts/self.h)
        if self.npts%self.h != 0:
            self.aligned = False
            self.aaa = np.vstack(tuple(self.aa for x in range(self.nhops + 1)))
        else:
            self.aligned = True
            self.aaa = np.vstack(tuple(self.aa for x in range(self.nhops)))

    def __call__(self, signal):
        a,g,err,h = self.lpc_fit(signal, self.order, self.h)
        return self.lpc_synth(a,g,err,h)
     


    def lpc_fit(self,signal):
        """
        Returns

        """
        signal = [float(x) for x in signal]
        signal = np.array(signal)
        if np.array(signal).ndim > 1:
            raise ValueError("Array of rank > 1 not supported")

        if self.order > signal.size:
            raise ValueError("Input signal must have a length >= LPC order")

        if self.order > 0:
            x = [i for i in signal]    
            p = self.order
            h = self.h
            w = 2*h
            npts = self.npts
            nhops = self.nhops
        
            edges = [0 for i in range(int((w-h)/2))]

            x = edges + x
            x = x + edges

            if self.aligned == False:
                a = np.zeros((nhops+1, p + 1))
                g = np.zeros(nhops+1)
                err = np.zeros(npts)
 
            else:
                a = np.zeros((nhops, p + 1))
                g = np.zeros(nhops)
                err = np.zeros(npts)

            for hop in range(nhops):
                xx = x[(hop)*h + 0: hop*h + w]
                aa = self.aa
                rs = scipy.signal.lfilter(aa, 1, xx[int((w-h)/2): int((w-h)/2) + h], axis = 0)
                G = math.sqrt(np.mean(rs**2))
                a[hop,:] = aa
                g[hop] = G
                err[(hop)*h + 0: hop*h + h] = rs/G
        
        if self.aligned == False:
            xx = signal[(self.nhops+1)*h:]
            rs = scipy.signal.lfilter(self.aa,1,xx)
            G = math.sqrt(np.mean(rs**2))
            a[hop + 1,:] = aa
            g[hop + 1] = G
            err[(hop+1)*h:] = rs/G
       
        self.a = a
        self.gains = g
        self.err = err        
  
        return a,g,err,self.npts,self.h
           
    def lpc_synth(self,a,g,err,npts,frame_width):
         
        e = err
        npts = npts
        d = np.zeros(npts)
        h = frame_width
        nhops = math.floor(npts/h)
        
        for hop in range(nhops):
            hbase = hop*h
            
            oldbit = d[hbase: hbase + h]
            aa = a[hop,:]
            
            G = g[hop]
           
            newbit = G*scipy.signal.lfilter(np.array([1.0]),aa,e[hbase : hbase + h])

            d[hbase : hbase + h] = newbit
        
        if self.aligned == False:
            hbase = (hop+1)*h
            oldbit = d[hbase:]
            #aa = a[hop+1,:]
            
            G = g[hop+1]
            
            newbit = G*scipy.signal.lfilter(np.array([1.0]),aa,e[hbase:])

            d[hbase:] = newbit

        self.d = d

        return d        

    def get_amp(self, err,frame):
        if self.aligned == True:
            self.amp = np.mean(err[::frame])
        else:
            self.amp = [np.mean(err[:self.nhops*self.h:self.h]),err[self.nhops*self.h]]
        if isinstance(self.amp,list):

            if self.amp[0] < 0:
                self.logger.debug("ERROR AMPLITUDE == NEGATIVE")
            else:
                self.logger.debug("ERROR AMPLITUDE == POSITIVE")

        else:
            if self.amp < 0: 
                self.logger.debug("ERROR AMPLITUDE == NEGATIVE")
            else:
                self.logger.debug("ERROR AMPLITUDE == POSITIVE")
        
        
        return self.amp
            
    def get_fits(self,err):
        h = self.h
        npts = self.npts
        nhops = math.floor(self.npts/self.h)
        fits = []
        if self.aligned == True:
            for i in range(nhops):
                p = np.polyfit(np.array(list(range((i*h)+2,(i+1)*h))),err[(i*h)+2:(i+1)*h],1)
                fits.append(p)

        else:
            for i in range(nhops):
                p = np.polyfit(np.array(list(range((i*h)+2,(i+1)*h))),err[(i*h)+2:(i+1)*h],1)
                fits.append(p)

            p = np.polyfit(np.array(list(range(nhops*h + 2,npts))),err[nhops*h + 2:],1)
            fits.append(p)

        self.fits = [x for sublist in fits for x in sublist]
        return self.fits


    def get_gains(self):
       return self.gains

    @classmethod
    def recon_err(cls,npts,frame_width,amp,fits):
    
        recon_error = np.zeros(npts)
        h = frame_width
        nhops = math.floor(npts/h)

        if npts%h == 0:       
            for i in range(nhops):
                recon_error[(i*h)+2:(i+1)*h] = np.polyval(np.array(fits[i]),np.array(list(range((i*h)+2,(i+1)*h))))        
            recon_error[::h] = amp
            recon_error[1::h] = -amp

        else:
            for i in range(nhops):
                recon_error[(i*h)+2:(i+1)*h] = np.polyval(np.array(fits[i]),np.array(list(range((i*h)+2,(i+1)*h))))
            recon_error[:npts-h:h] = amp[0]
            recon_error[1:npts-h:h] = -amp[0]

            recon_error[nhops*h] = amp[1]
            recon_error[nhops*h + 1] = -amp[1]
            recon_error[nhops*h + 2:] = np.polyval(np.array(fits[i+1]),np.array(list(range(nhops*h + 2,npts))))
        cls.r_err = recon_error

        return recon_error
    
    
    def pack_residual(self,f):
        packed = []
        packed.append(self.gains)
        packed.append(self.fits)
       
        
        f.write(struct.pack('@I',self.npts))
        f.write(struct.pack('@I',self.h))

        if self.aligned == True:
            
            f.write(struct.pack('@d',self.amp))
        else:
            for item in self.amp:
                f.write(struct.pack('@d',item))

        for item in packed:
            for x in item:
               f.write(struct.pack('@d', x))
        
        
        
        return f
    
    @classmethod    
    def unpack_residual(cls,f):
        
        f.seek(0)
 
        npts = struct.unpack('@I',f.read(struct.calcsize('I')))[0]
        frame_size = struct.unpack('@I',f.read(struct.calcsize('I')))[0]
        nhops = math.floor(npts/frame_size)
        
        if npts%frame_size == 0:
            amp = struct.unpack('@d',f.read(struct.calcsize('d')))[0]
            gains = []
            for i in range(nhops):
               gains.append(struct.unpack('@d',f.read(struct.calcsize('d')))[0])

        else:
            amp = []
            for i in range(2):
                amp.append(struct.unpack('@d',f.read(struct.calcsize('d')))[0])
            gains = []
            for i in range(nhops + 1):
                gains.append(struct.unpack('@d',f.read(struct.calcsize('d')))[0])
        
        
        fits = []
        '''
        while True:
        
            try:
                fits.append(struct.unpack('@d',f.read(struct.calcsize('d')))[0])

            except struct.error:
                print(fits
                fits = np.reshape(fits,(len(gains),2))
                fits = fits.tolist()                     
                break
        '''
        for i in range(2):
            fits.append(struct.unpack('@d',f.read(struct.calcsize('d')))[0])

        fits = np.reshape(fits,(1,2))
        fits = fits.tolist()
        
        return npts, frame_size, amp, gains, fits, f
       
                

if __name__ == "__main__":
    data = loadmat('r10.mat')
    data = data['residual'][0]
    ts = data[:1400]
    #ts = list(range(100))
    #a,g,err,h = lpc_fit(ts, 2, 350)
    #print(a)
    lpc = LPC(2,300,1400)
    a,g,err,npts,h = lpc.lpc_fit(ts)
    lpc.lpc_synth(a,g,err,npts,h)
    print(lpc.npts)
    print(lpc.get_amp(err,lpc.h))
    lpc.get_gains()
    lpc.get_fits(err)
    lpc.pack_residual()
    npts,frame_width,amp,gains,fits = lpc.unpack_residual()
    recon_err = lpc.recon_err(npts,frame_width,amp,fits)
    #print((recon_err - lpc.err).tolist())
    r_err = lpc.lpc_synth(lpc.aaa,gains,LPC.r_err,npts,frame_width)
    print(ts-r_err)
    
