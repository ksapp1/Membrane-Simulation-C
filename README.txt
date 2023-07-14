VESICLE_SIM

INTRODUCTION:

vesicle_sim is a C program that performs numerical simulations of vesicles, either approximated as a plane or a sphere, coupled to the diffusion curvature dependent lipids.

REQUIREMENTS:

GSL libraries are required for random number generation and spherical harmonics.

INSTALLATION:

autoreconf --install
configure
make

If the automake files are changed you may need to run: autoreconf

program is built in "./src/", to change build location update the automake files: configure.ac and Makefile.am

USAGE:

./src/vesicle_sim PARAM_FILE

See example parameters file ("./src/parametrs.inp") for syntax in the parameter input file


