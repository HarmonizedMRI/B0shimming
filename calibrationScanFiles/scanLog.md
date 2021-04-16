

## UM MR750 3T (inside scanner)


### Calibration

8ch head array  
FBIRN phantom  
toppev3  
Exam 12775, Series 4

Auto prescan: x = -1, y = 2, z = 1, cf = 127788922

'x', 'y', 'z', 'z2', 'xy', 'zx', 'x2y2', 'zy'

Scans:  
system shims  
z2 500     (setNavShimCurrent z2 500)
z2 -500  
xy 500  
xy -500  
zx 500  
zx -500  
x2y2 500  
x2y2 -500  
zy 500  
zy -500  
x 9     (manual prescan)  
x -11  
y 12  
y -8  
z 11  
z -9  
system shims


## Move scanArchive files from scanner to quickstep:
```
$ cd /export/home1/raw_temp/jfnielse/
$ rsync -avz 'sdc@vre:/data/arc/Closed/Exam12775/Series4/*.h5' .
$ rsync -avz . jfnielse@quickstep:/mnt/storage/jfnielse/B0shimming/UM-10-Dec-2020/


