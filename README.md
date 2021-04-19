# JPetalo

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> JPetalo

It is authored by jjgomezcadenas.

To (locally) reproduce this project, do the following:

-1. Install Julia from https://julialang.org/ 
0.  Download this code base. The raw data is not included in the
   git-history, see below for instructions on how to dowload the example file.
   
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages.

The interesting stuff, so far, is in notebook recolor.jl (reconstruction LOR). recolor loads a databa with the positions of SiPMs in a full-body pet (fbpet) and an example file with the charge deposited in the SiPMs by
two gamma events. It then clasifies the SiPMs in each event (as being in the same side or in the opposite side
of the SiPM with more charge, thus dividing the detector in two), clusters the SiPMs (after a cut in charge) and
computes the line of response. 

2. To run the recolor.jl Pluto notebook you need Pluto:
   ```
   julia> import Pluto
   julia> Pluto.run()
   ```

Notice that recolor.jl is a reactive notebook. All cells will run before you get control.
This can take a little, while everything compiles. Now is your change to prepare that coffee.

3. Dr. Watson defines a few useful goodies, including datadir() that points to the directory where you
can keep your files. recolor.jl assumes that the PETALO SiPM database is in
```
datadir("db/petalodb.csv")
```
And also will search for a file
```
datadir("fbpet/full_body_phantom_paper.19.h5.csv")
```

You can download both from the directory:
```
http://next.ific.uv.es/icgallery/public/PETALO
```

You can click in the files and save them in their respective directories (db and fbpet).

recolor is a live notebook. Try a few events and modifying the energy cut. Enjoy!



