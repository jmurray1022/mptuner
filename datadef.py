import csv
import wx
from support import *
import numpy as np
import matplotlib.numerix.ma as ma

class DataObj:
    def __init__(self, fname):
        self.fname=fname
        self.ind=array([])
        self.data=array([])
        self.colnames=array([])
        self.subset=array([])
        self.cols=array([])
        self.colsettings=array([])
        self.factors=[]
        self.settings=[]
        self.groups=[]
        self.groupnames=[]
        self.group_cols=[]
        self.groupon=[]
        
        self.read_data()
        self.init_masks()
    
    def read_data(self):
        '''Read the data in self.fname. Self.colnames contains column names, self.data has the data as floats.
        '''
        
        file=open(self.fname)
        line=file.readline()
        if line[-1]=='\n':
            line=line.rstrip()
        if line[0]=='"':
            line=line.replace('"','')
        cn=line.rsplit(",")
        self.indname = cn[0]
        self.colnames=np.array(cn[1:len(cn)])
        nc=len(self.colnames)
        self.colsettings=np.vstack([name.rsplit(".")] for name in self.colnames)
        #self.settings

        def missing_data(x):
            if x == 'NaN' or x=='NA':
                return np.nan
            else:
                return float(x)
        conv = dict( (k, missing_data) for k in range(nc+1) )
        df = np.loadtxt(self.fname, delimiter = ',', converters = conv, skiprows=1)
        self.ind=df[:,0]
        self.data=df[:,1:len(df[0,])]
        del df
        #self.data=ma.masked_invalid(df[:,1:len(df[0,])])
        
        self.factors = map(str,range(len(self.colsettings[0])))
        self.settings = dict((self.factors[i], sorted(list(set(self.colsettings[:,i])), ncmp)) for i in range(len(self.factors)))
        
    def init_masks(self):
        self.resmask = ma.make_mask_none(self.data.shape)
        self.colmask = ma.make_mask_none(self.data.shape)
    
    def do_seq_polish(self, frame, cols=array([]), cs=1.5, ce=3.0, eps=0.01, maxiter=100, repolish=True):
        '''Do the sequential median polish. Build two masks, one for coleff and one for reseff.
        Required: 
        frame, int representing frame size
        cols, numpy array with column indices
        '''
        shape=self.data.shape
        nrow=shape[0]
        ncol=shape[1]

        if cols==array([]):
            cols=arange(0,ncol)
            
        numframe=nrow/frame
        r=nrow%frame
        
        print r; print numframe
        
        #The data we skip due to frame size is all masked
        self.colmask[0:r,:]=ones(self.colmask[0:r,:].shape,dtype=bool)
        self.resmask[0:r,:]=ones(self.resmask[0:r,:].shape,dtype=bool)
        
        dlg = wx.ProgressDialog("Median polish", "Performing median polish algorithm", 
            maximum = numframe, style = wx.PD_APP_MODAL | wx.PD_SMOOTH | wx.PD_AUTO_HIDE)
  
        for i in range(0,numframe):
            low=r+i*frame
            high=r+(i+1)*frame
            polish=medpolish(self.data[low:high,cols], eps=eps, maxiter=maxiter)
            self.colmask[low:high,cols]=colclean(self.data[low:high,cols], polish['col'], cs)
            self.resmask[low:high,cols]=resclean(self.data[low:high,cols], polish['res'], ce)
            dlg.Update(i+1)
            
        dlg.Destroy()
    def outdata(self):
        mask = ma.mask_or(self.colmask, self.resmask)
        #mask[0,1] = 1
        #mask[0,2] = 1
        #mask[0,4] = 1
        #clean = ma.masked_array(data=self.data, mask = ma.mask_or(self.colmask, self.resmask), fill_value = "NA")
        clean = ma.masked_array(data=self.data, mask = mask, fill_value = "NA")
        #np.savetxt('out.csv', clean.filled(np.nan), delimiter = ",")
        ldata = clean.filled(np.nan).tolist()
        for i, row in enumerate(ldata):
            for j, col in enumerate(row):
                if np.isnan(ldata[i][j]): ldata[i][j] = "NA"
        fulldata = vstack((hstack((array(self.indname),self.colnames)), hstack((transpose(np.atleast_2d(self.ind)), ldata))))
        #writer = csv.writer(open("out2.csv", "wb"))
        #writer.writerows(fulldata)
        return fulldata
#data=DataObj("data.csv")
#data.read_data()
#data.init_masks()
#data.do_seq_polish(8)
#data.polish()