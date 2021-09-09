Model GPS receiver. 
Uses trilateration to compute location given satellite data. 
The (Java) vehicle uses `vehicle.class` and `angles.class`. 
To compile satellite and receiver programs, run
```
gcc gps.c receiver.c matops.c cholSolve.c getSatellite.c -o receiver
gcc gps.c satellite.c matops.c cholSolve.c getSatellite.c -o satellite
```
(`gps.c` contains helper functions for both satellite and receiver and `matops.c` and `qrSolve.c` contain linear algebra operations.)
To run, pipe a trip file, e.g. `bm.dat`, into `vehicle` then to the satellite then to the receiver:
```
cat bm.dat | java vehicle | ./satellite | ./receiver
```
This should approximate the location data given by the vehicle as in the command line 
```
cat bm.dat | java vehicle 
```

TODO: implement Cholesky decomposition instead of QR to exploit symmetry of Jacobian in receiver. 
