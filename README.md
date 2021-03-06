# Stochastic Climate

A reproduction of the model in: [https://arxiv.org/pdf/2104.06216.pdf](https://arxiv.org/pdf/2104.06216.pdf)

-----

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> Stochastic Climate

It is authored by Mark Baum <markmbaum@protonmail.com> and Kaitlyn Loftus <kloftus@g.harvard.edu>.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
