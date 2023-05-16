import sys,os
dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.insert (0, dir_path)
os.chdir(dir_path)
from app import app as application
