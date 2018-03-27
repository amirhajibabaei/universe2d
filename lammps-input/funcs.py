


def configure():
   import os
   import glob
   configfile = 'configure.txt'
   file1 = 'therm_*.txt'
   file2 = 'snap_*.txt'
   os.system('rm -f '+configfile)
   config = None
   if   os.path.isfile('restart1'): 
      config = 'restart1'
   elif os.path.isfile('restart2'): 
      config = 'restart2'
   elif len(glob.glob('restart_*'))>0:
      config = 'restart_*'
   if not config:
      os.system('cat ini.box > '+configfile)
      #os.system('rm -f '+file1+' '+file2)
      return 1
   with open(configfile,'w') as f:
      f.write('read_restart '+config+'\n')
   index1 = max(int(v.replace('therm_','').replace('.txt','')) for v in glob.glob(file1))
   index2 = max(int(v.replace('snap_','').replace('.txt','')) for v in glob.glob(file2))
   index  = max(index1, index2)
   return index + 1


