gfortran -coverage -fprofile-arcs -ftest-coverage -O0 ./cosmic/src/benchmarkevolv2.f ./cosmic/src/comenv.f ./cosmic/src/corerd.f ./cosmic/src/deltat.f ./cosmic/src/dgcore.f ./cosmic/src/evolv2.f ./cosmic/src/gntage.f ./cosmic/src/hrdiag.f ./cosmic/src/instar.f ./cosmic/src/kick.f ./cosmic/src/mix.f ./cosmic/src/mlwind.f ./cosmic/src/mrenv.f ./cosmic/src/ran3.f ./cosmic/src/rl.f ./cosmic/src/star.f ./cosmic/src/zcnsts.f ./cosmic/src/zfuncs.f ./cosmic/src/concatkstars.f ./cosmic/src/bpp_array.f ./cosmic/src/checkstate.f ./cosmic/src/eddington.f ./cosmic/src/pisn.f ./cosmic/src/fallback.f -o benchmarkevolv2.exe -I cosmic/src -Wl,-rpath,${CONDA_PREFIX}/lib
./benchmarkevolv2.exe
