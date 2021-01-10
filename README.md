Model GPS receiver. 
Uses trilateration to compute location given satellite data. 
The (professor given, Java) vehicle uses `vehicle.class` and `angles.class`. 
To run, compile the `C` code as
```
gcc satellite.c -o sat
gcc receiver.c -o rec
```
Then, pipe a trip file, e.g. `bm.dat`, into `vehicle` then to the satellite then to the receiver:
```
cat bm.dat | java vehicle | sat | rec
```
This should approximate the same location data given by the vehicle as in the command line 
```
cat bm.dat | java vehicle 
```
