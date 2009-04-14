from pylab import *
import numpy as np
import matplotlib
#import matplotlib.dates as mdates
import matplotlib.mlab as mlab
import matplotlib.numerix.ma as ma
import csv

#class Group

def ncmp(a, b):
    try:
        x=float(a)
        y=float(b)
    except ValueError:
        x=a
        y=b
    if x>y:
       return 1
    elif x==y:
       return 0
    else: # x<y
       return -1


###########################################################
# Begin quantile function
###########################################################

def quantile(a, p):
    """Returns the pth quantile of a (a numpy array). Interpolation is 
    type 7 described in the R documentation, the default for R's quantile()
    """
    x=a[where(~isnan(a))]
    if not x.any():
        return np.nan
    else:
        n=len(x)
        x.sort()
        k=p*(n-1)+1
        i=int(k) #truncates k toward zero, gives us an index
        g=k-i #the fractional part of k, used to weight the 2 closest order statistics
        if i<n:
            result=(1-g)*x[i-1]+g*x[i]
        else:
            result=x[i-1]
    return result

###########################################################
# End quantile function
###########################################################

###########################################################
# Begin iqr (interquartile range) function
###########################################################

def iqr(x):
    '''Returns the interquartile range of a numpy array x.
    '''
    result=quantile(x,0.75)-quantile(x,0.25)
    return result

###########################################################
# End iqr (interquartile range) function
###########################################################

###########################################################
# Begin medpolish (median polish) function
###########################################################

def medpolish(z,eps=0.01,maxiter=10):
    '''Perform Tukey's median polish on z with precision eps.
    The algorithm will continue until convergence or until maxiter'''
    #Initializing variables.
    x=z.copy()
    nr=x.shape[0]
    nc=x.shape[1]
    t=0
    delta=0
    r=np.zeros(nr)
    c=np.zeros(nc)
    rdelta=r
    cdelta=c
    oldsum=np.sum(x)
    newsum=0
    converged=0 #0==false
    result="Error"

    #begin polish loop
    for i in range(1,maxiter):

        #subtract the row medians from the rows of the matrix
        #and add them to r
        np.median(x,1)
        x=x-rdelta[:,newaxis]
        r=r+rdelta

        #subtract the median of r from r and add it to the overall effect
        delta=np.median(r)
        r=r-delta
        t=t+delta

        #subtract the column medians from the columns of the matrix
        #and add them to c
        cdelta = np.median(x,0)
        x=x-cdelta
        c=c+cdelta

        #subtract the median of c from c and add it to the overall effect
        delta=np.median(c)
        c=c-delta
        t=t+delta

        #sum of abs value of the elements
        newsum=np.sum(np.abs(x[where(~isnan(x))]))

        #the polish has converged if the abs sum of the residual matrix
        #is arbitrarily small (or zero)
        converged=abs(newsum-oldsum)<eps*newsum or newsum==0

        if converged:
            break

        oldsum=newsum
    #End polish loop

    if converged:
        #make our results
        result={"iter":i-1, "overall":t, "row":r, "col":c, "res":x}
    else:
        print "median polish failed to converge. increase maxiter?"
        #result={"iter":i-1, "overall":t, "row":r, "col":c, "res":x}
    return result

##########################################################################
# End medpolish function
##########################################################################

##########################################################################
# Begin colclean (column effect clean) function
##########################################################################

def colclean(x, col_eff, cs):
    iqrange=iqr(col_eff)
    low=quantile(col_eff,0.25)-cs*iqrange
    high=quantile(col_eff,0.75)+cs*iqrange
    colkill=where((col_eff<low)|(col_eff>high))
    mask=ma.make_mask_none(x.shape)
    mask[:,colkill]=ones(mask[:,colkill].shape,dtype=bool)
    return mask
    
def rcolclean(x, col_eff, cs):
    iqrange=iqr(col_eff)
    low=quantile(col_eff,0.25)-cs*iqrange
    high=quantile(col_eff,0.75)+cs*iqrange
    colkill=where((col_eff<low)|(col_eff>high))
    mask=ma.make_mask_none(x.shape[1])
    mask[colkill]=ones(mask[colkill].shape, dtype=bool)
    return mask
    
##########################################################################
# End colclean (column effect clean) function
##########################################################################

##########################################################################
# Begin resclean (residuals clean) function
##########################################################################
def resclean(x, res, ce):
    y=x.ravel()
    iqrange=iqr(y)
    low=quantile(y,0.25)-ce*iqrange
    high=quantile(y,0.75)+ce*iqrange
    mask=ma.masked_outside(x,low,high).mask
    #if mask.dtype == np.dtype('bool'):
    #    mask = zeros(x.shape,dtype=bool)
    return mask
    
def rresclean(x, lag, res, ce):
    y=x.ravel()
    iqrange=iqr(y)
    low=quantile(y,0.25)-ce*iqrange
    high=quantile(y,0.75)+ce*iqrange
    mask=ma.masked_outside(x,low,high).mask
    #if mask.dtype == np.dtype('bool'):
    #    mask = zeros(x.shape,dtype=bool)
    return mask[lag,:]
    
#http://automatthias.wordpress.com/2007/04/28/cartesian-product-of-multiple-sets/

def cprod(lists, previous_elements = []):
    if len(lists) == 1:
        for elem in lists[0]:
            yield previous_elements + [elem, ]
    else:
        for elem in lists[0]:
            for x in cprod(lists[1:], previous_elements + [elem, ]):
                yield x