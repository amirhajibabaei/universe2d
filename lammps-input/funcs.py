


def configure():
   import os
   import glob
   import time
   configfile = 'configure.txt'
   file1 = 'md_vectors_*.txt'
   file2 = 'md_snap_*.txt'
   os.system('rm -f '+configfile)
   config = None
   if   os.path.isfile('md_restart1'): 
      config = 'md_restart1'
   elif os.path.isfile('md_restart2'): 
      config = 'md_restart2'
   elif len(glob.glob('md_restart_*'))>0:
      config = 'md_restart_*'
   if not config:
      os.system('cat ini.box > '+configfile)
      #os.system('rm -f '+file1+' '+file2)
      return 1
   with open(configfile,'w') as f:
      f.write('read_restart '+config+'\n')
   index1 = max(int(v.replace('md_vectors_','').replace('.txt','')) for v in glob.glob(file1))
   index2 = max(int(v.replace('md_snap_','').replace('.txt','')) for v in glob.glob(file2))
   index  = max(index1, index2)
   time.sleep(1)
   return index + 1


