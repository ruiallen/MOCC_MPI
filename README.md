# MOCC_MPI
Full quantum dynamic versions of QMOCC, with MPI integrated. Modified based on Stancil's and CY-lin's version. 
QMOCC_MPI is parallelized on the number of partial waves used for expansion, evenly distributed.
Requires 
  1, input.capture file, which includes calculation parameters (like dr and Rmax) and all state information
  2, pot.diabatic which is the coupling matrix with potentials and coupling information (Antisymmetric Matrix)

outputs:
  1, all .txt files contain information about scattering amplitude.
  2, rmat and imat are real part and imaginary parts of the scattering S-matrix respectively. 

To compile it, use the Makefile provided. Tested with intel's mpiifort with proper flags. 
to run, simply do srun ./mocc. The node number and the number of core information have already been taken care of inside the code.
