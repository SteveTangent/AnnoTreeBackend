import sys,os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert (0, dir_path)
os.chdir(dir_path)

#f = open('/datadrive/AnnoTree/AnnoTreeBackend/wsgilog.txt', 'a')
#timestr = time.strftime("%Y.%m.%d-%H:%M:%S")
#f.write(timestr+'\n')
#f.write(sys.version)
#f.write('running\n')
#f.close()

from app import app as application
