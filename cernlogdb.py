"""
Wrapper around the  CERN logging database command line tool

Usage:
  mkcmd:  Produce a string for the CERN logging database command line tool
  getdb:  Query and parse the results
  open:   Parse and return the output of a query from a file name
  load:   Parse and return the output of a query from a file object
  pprint: Describe the output of getdb, open, load
"""


import os
import time
from glob import glob

import numpy as _np

import rdmdate as _date

import time as _t


def mkcmd(vs,t1,t2,step=None,scale=None,exe='cern-logdb',conf='ldb.conf',filename=None,format=None):
  """Produce a string for the CERN logging database command line tool.

  Usage getdata("SPS.BCTDC.31832:INT","2010-06-10 00:00:00","2010-06-10 23:59:59",step='20 SECOND')

  Arguments:
    vs:      name of the variable in the database
    t1,t2:  start and end time as string in the %Y-%m-%d %H:%M:%S.SSS format
            or unix time
    step:   For multiple file request '<n> <size>'
    scale:  For scaling algoritm '<n> <size> <alg>'
    format: For output file format in CSV, XLS, TSV, MATHEMATICA
    filename: String for output name

  where:
    <n> is an integer number
    <size> is one of SECOND, MINUTE, HOUR, DAY, WEEK, MONTH, YEAR
    <alg> one of AVG, MIN, MAX, REPEAT, INTERPOLATE, SUM, COUNT
    <date>
  """
  if type(t1) in [float,int]:
    if t1<0:
      t1+=_t.time()
    t1=_date.dumpdate(t1)
  if type(t2) in [float,int]:
    if t2<0:
      t2+=_t.time()
    t2=_date.dumpdate(t2)
  cmd='%s -vs "%s" -t1 "%s" -t2 "%s"' %(exe,vs,t1,t2)
  if conf:
    cmd+=' -C %s' % conf
  if scale:
    n,size,alg=scale.split()
    cmd+=' -sa "%s" -ss "%s" -si "%s"' %(alg,n,size)
  if step:
    ni,it=step.split()
    cmd+=' -IT "%s" -NI "%s"' %(it,ni)
  if filename:
    cmd+=' -N "%s"' % filename
  if format:
    cmd+=' -F "%s"' % filename
  return cmd

_interval={'day': 86400,
       'hour': 3600,
       'minute': 60,
       'month': 2592000,
       'second': 1,
       'week': 604800,
       'year': 31536000}


def generate_filenames(mask,t1,t2,step):
  """ Generate filenames as done by the CERN DB command line tool

  Arguments:
    mask string
    step:   For multiple file request '<n> <size>'
    t1,t2:  start and end time as string in the %Y-%m-%d %H:%M:%S.SSS format
            or unix time

  where:
    <n> is an integer number
    <size> is one of SECOND, MINUTE, HOUR, DAY, WEEK, MONTH, YEAR
  """
  if type(t1) in [float,int]:
    t1=_date.dumpdate(t1)
  if type(t2) in [float,int]:
    t2=_date.dumpdate(t2)
  ni,it=step.split()
  sstep=int(ni)*_interval[it]



def _floatarray(i):
  return _np.array(i,dtype=float)


def load(fh,sep=',',ttype=float,vtype=_floatarray,debug=False):
  """Parse the output of the CERN logging database command line tool

  Usage:
    fh:    file handler
    sep:   separator ',' or '\\t'
    ttype: function to convert the timestamp
    vtype: function to covernt the record

  Returns:
    A dictionary for which for each variable found data is stored in a tuple
    of timestamps list and record list. Data is accesible by the variable name
    and a numeric index.
    In addition the keys:
      'log'      contains the log messages contained in the file
      'datavars' contains the a list of variable names
  """
  data={}
  dataon=False
  header=True
  log=[]
  datavars=[]
  for l in fh:
    if l.startswith('VARIABLE'):
      vname=l.split()[1]
      #print 'Found var %s' % vname
      if vname in data:
        t,v=data[vname]
      else:
        t,v=[],[]
        data[vname]=(t,v)
        datavars.append(vname)
      dataon=False
      header=False
    elif l.startswith('Timestamp'):
      dataon=True
    elif l=='\n':
      dataon=False
    elif dataon:
      ll=l.strip().split(sep)
      t.append(ttype(ll[0]))
      v.append(vtype(ll[1:]))
    elif header:
      if debug is True:
        print l
      log.append(l)
  if len(log)>0:
    data['log']=log
  data['datavars']=datavars
  for ii,nn in enumerate(datavars):
    data[ii]=data[nn]
  return data

def combine_data(data,vtype=float,ttype=int):
  """Combine vector data"""
  for ik,k in enumerate(data['datavars']):
    t,v=data[k]
    data[k]=[_np.array(t,dtype=ttype),_np.array(v,dtype=vtype)]
    data[ik]=data[k]
  return data

def open(fn,sep=','):
  """Load output of the CERN measurement database query from a filename

  Usage:  open("test.tsv")
    fn: filename
    sep: separator type

  Separator is inferred from the file extension as well
  """
  fh=file(fn)
  if fn.endswith('tsv'):
    sep='\t'
  data=load(fh,sep=sep)
  fh.close()
  data['filename']=fn
  return data

def openglob(mask,sep=','):
  """Open a set of filenames"""
  fnames=iglob(mask)
  data=load(_icat(fnames),sep=sep)
  data['filename']=fnames
  return data

def _icat(fnames):
  for fn in fnames:
    for l in file(fn):
      yield l


def dbget(vs,t1,t2,step=None,scale=None,exe='cern-logdb',conf='ldb.conf',
         ttype=str,vtype=list):
  """Query the CERN measurement database and return data
  Usage dbget("SPS.BCTDC.31832:INT","2010-06-10 00:00:00","2010-06-10 23:59:59",step='20 SECOND')

  Arguments:
    vs:      name of the variable in the database
    t1,t2:  start and end time as string in the %Y-%m-%d %H:%M:%S.SSS format
            or unix time
    step:   For multiple file request '<n> <size>'
    scale:  For scaling algoritm '<n> <size> <alg>'

  where:
    <n> is an integer number
    <size> is one of SECOND, MINUTE, HOUR, DAY, WEEK, MONTH, YEAR
    <alg> one of AVG, MIN, MAX, REPEAT, INTERPOLATE, SUM, COUNT
    <date>
  """
  cmd=mkcmd(vs,t1,t2,step,scale,exe,conf)
  fh=os.popen(cmd)
  data=load(fh,ttype=ttype,vtype=vtype)
  fh.close()
  return data


def pprint(data):
  """Pretty print data, last dimension is from the first record"""
  print "CERN DB data:"
  for k in data['datavars']:
    t,v=data[k]
    recl=set()
    for ii,vv in enumerate(v):
      recl.add(len(vv))
    recl=' or '.join([str(i) for i in recl])
    print "  ['%s'][1] => v[%d][%s]" %( k,len(t),recl)


def dbget_repeat(vs,t1,t2,step,sa=None,sf=None,exe='cern-ldb'):
  if type(t1) is str:
    t1=round(_date.parsedate_myl(t1),0)
  if type(t2) is str:
    t2=round(_date.parsedate_myl(t2),0)
  data_final={}
  data_final['log']=[]
  data_final['datavars']=vs.split(',')
  for vname in data_final['datavars']:
    data_final[vname]=([],[])
  for t in range(int(t1),int(t2),step):
    nt1=_date.dumpdate(t)
    nt2=_date.dumpdate(t+step-0.001)
    print 'calling dbget from %s to %s' % (nt1,nt2)
    data=dbget(vs,nt1,nt2,sa=sa,sf=sf,exe=exe)
    for vname in data['datavars']:
      t,v=data[vname]
      nt=data_final[vname][0].extend(t)
      nv=data_final[vname][0].extend(v)
    data_final['log'].extend(data['log'])
  return data_final

def merge_out(fnames):
  data_final={}
  data_final['log']=[]
  data_final['datavars']=vs.split(',')
  for vname in data_final['datavars']:
    data_final[vname]=([],[])
  for fn in fnames:
    data=load(fn)
    for vname in data['datavars']:
      t,v=data[vname]
      nt=data_final[vname][0].extend(t)
      nv=data_final[vname][1].extend(v)
    data_final['log'].extend(data['log'])
  return data_final



if __name__=='__main__':
  print mkcmd('CMS:INTLUMI_DELIVERED_STABLEBEAMS',
        t1='2010-09-01 00:00:00.000',
        t2='2010-09-30 23:59:59.999',scale='1 DAY MAX')


