import sys
import argparse

from collections import OrderedDict

from mip_pipe_common import *

def main():
    opts = argparse.ArgumentParser( description='clip off MID from MIP reads and place in a separate file' )

    opts.add_argument('--inFq', default=None, dest='inFq')
    opts.add_argument('--midLen', default=10,type=int, dest='midLen')
    opts.add_argument('--outMidList', default=None, dest='outMidList')
    opts.add_argument('--outFq', default=None, dest='outFq')

    o = opts.parse_args()

    filInFq = zopen(o.inFq,'r')

    ML = o.midLen

    if o.outMidList is not None and o.outFq is not None:
        filOutMidList = zopen(o.outMidList,'w')
        filOutFq = zopen(o.outFq,'w')

        l=filInFq.readline()
        while len(l)>0:
            filOutFq.write(l)
            filOutMidList.write(l)

            l=filInFq.readline().rstrip() 
            
            filOutMidList.write('%s\n'%l[-ML:])
            filOutFq.write('%s\n'%l[:-ML])

            l=filInFq.readline()

            filOutFq.write(l)
            filOutMidList.write(l)

            l=filInFq.readline().rstrip() 

            filOutMidList.write('%s\n'%l[-ML:])
            filOutFq.write('%s\n'%l[:-ML])
            
            l=filInFq.readline()
        
        filOutMidList.flush()
        filOutMidList.close()

        filOutFq.flush()
        filOutFq.close()

    elif o.outMidList is None:
        filOutFq = zopen(o.outFq,'w')

        l=filInFq.readline()
        while len(l)>0:
            filOutFq.write(l)
            l=filInFq.readline().rstrip()
            
            filOutFq.write('%s\n'%l[:-ML])
            l=filInFq.readline()

            filOutFq.write(l)
            l=filInFq.readline().rstrip()

            filOutFq.write('%s\n'%l[:-ML])            
            l=filInFq.readline()
        
        filOutFq.flush()
        filOutFq.close()

    elif o.outFq is None:
        filOutMidList = zopen(o.outMidList,'w')

        l=filInFq.readline()
        while len(l)>0:
            filOutMidList.write(l)

            l=filInFq.readline().rstrip()
            
            filOutMidList.write('%s\n'%l[-ML:])

            l=filInFq.readline()

            filOutMidList.write(l)

            l=filInFq.readline().rstrip()

            filOutMidList.write('%s\n'%l[-ML:])
            
            l=filInFq.readline()
        
        filOutMidList.flush()
        filOutMidList.close()

if __name__ == '__main__':                
    main()