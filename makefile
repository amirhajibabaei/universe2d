mc_exe: universe2d.f08 universe2d_mpi.f08 feeder.f08 metropolis_mpi.f08
	mpifort -o mc_exe -O3 -flto -Wall -Wextra universe2d.f08 universe2d_mpi.f08 feeder.f08 metropolis_mpi.f08
