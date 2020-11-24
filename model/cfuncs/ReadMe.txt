////////////////////
// File structure //
////////////////////

MEX gateways to be called from MATLAB
includes/:          .cpp files
nlopt-2.4.2-dll64/: files for NLopt (unpacked from nlopt-2.4.2-dll64.zip)


//////////////
// GATEWAYS //
//////////////

For solving:      solve_model.cpp
For simulating:   simulate_model.cpp

Output (in MATLAB folder):
log_solve.txt
log_simualte.txt
log_assert.txt


//////////////
// CONTROLS //
//////////////

How much to print in logs? See macros in includes/logs.cpp
How many assertions to test? See macros in includes/assert.cpp