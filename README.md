Model GPS receiver. 
Uses trilateration to compute location given satellite data. 
Location and time data are given in the following format: 
<time> <latitude minutes> <latitude degrees> <latitude seconds> <northern/southern hemisphere> <longitude minutes> <longitude degrees> <longitude seconds> <eastern/western hemisphere> <elevation>
The (Java) vehicle uses `vehicle.class` and `angles.class`. 
To compile satellite and receiver programs, run
```
gcc gps.c receiver.c matops.c cholSolve.c getSatellite.c -o receiver
gcc gps.c satellite.c matops.c cholSolve.c getSatellite.c -o satellite
```
`gps.c` contains helper functions for both satellite and receiver and `matops.c` and `qrSolve.c` contain linear algebra operations. 
Linear algebra files can be found in nla repository. 
To run, pipe a trip file, e.g. `bm.dat`, into `vehicle` then to the satellite then to the receiver. 
The following uses a trip from Salt Lake City, UT to the point of zero longitude and latitude:
```
$ cat o.dat | java vehicle | ./satellite | ./receiver
102123.210000 40 45 55.000155 1 111 50 58.000095 -1 1371.999970
103220.630000 36 14 8.880108 1 99 25 18.220047 -1 1219.548082
104318.050000 31 42 22.769938 1 86 59 38.439948 -1 1067.111167
105415.470000 27 10 36.659983 1 74 33 58.659937 -1 914.660985
106512.890000 22 38 50.550094 1 62 8 18.880087 -1 762.217440
107610.310000 18 7 4.440011 1 49 42 39.109978 -1 609.770477
108707.730000 13 35 18.330040 1 37 16 59.330011 -1 457.326804
109805.150000 9 3 32.220014 1 24 51 19.550042 -1 304.877032
110902.570000 4 31 46.109995 1 12 25 39.770018 -1 152.441257
112000.000000 0 0 0.000032 1 0 0 0.000038 -1 0.000163
```
This should approximate the location data given by the vehicle as in the command line 
```
$ cat o.dat | java vehicle 
102123.21 40 45 55.0 1 111 50 58.0 -1 1372.0
103220.63 36 14 8.88 1 99 25 18.22 -1 1219.55
104318.05 31 42 22.77 1 86 59 38.44 -1 1067.11
105415.47 27 10 36.66 1 74 33 58.66 -1 914.66
106512.89 22 38 50.55 1 62 8 18.88 -1 762.22
107610.31 18 7 4.44 1 49 42 39.11 -1 609.77
108707.73 13 35 18.33 1 37 16 59.33 -1 457.33
109805.15 9 3 32.22 1 24 51 19.55 -1 304.88
110902.57 4 31 46.11 1 12 25 39.77 -1 152.44
112000.0 0 0 0.0 1 0 0 0.0 1 0.0
```
