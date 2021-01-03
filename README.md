Model GPS receiver. 
Uses trilateration to compute location given satellite data. 
The (professor given) vehicle uses `vehicle.class` and `angles.class`. 
To run, compile the `C` code as
```
gcc sat.c -o sat
gcc rec.c -o rec
```
**NOTE**: the executable names `satellite` and `receiver` are already taken by the professor given programs. 

Then, pipe a trip file, e.g. `bm.dat`, into `vehicle` then to your satellite then to your receiver:
```
cat bm.dat | java vehicle | sat | rec
```
This should approximate the same location data given by the vehicle as in the command line 
```
cat bm.dat | java vehicle 
```
