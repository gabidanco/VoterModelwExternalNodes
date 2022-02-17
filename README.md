This program is written in fortran 90 and simulates the voter model
dynamics with external varying nodes in regular networks.

The dynamics for the external nodes are such that the focal node augments
with 'lambda' its external influence aligned with its neighborhood, and
diminishes by 'lambda' their contrary external influence. Focal does not 
consider its own state, only the average opinion of their neighbors.

Any suggestions, comments or questions can be sent to:
gabidanco@gmail.com.

________________________________________________________________________

To run the program:

1) Compile the file 'voter_dyn_ext.f90' in your operating system.
2) Put the executable file and the input.in file in the same directory.
3) Run the executable file 'voter_dyn_ext.exe'

--------Input files-----------

---> "seed.in"
- this is a file with 12 integer numbers separated by blank spaces.
- this file is updated at the end of each run.

---> "input_de.in"
- n: total number of dynamic nodes
- deg: degree of the network will be 2*deg
- np: total number of external influence, to be kept constant, np=n0+n1
- tau: equilibrium time 
- nm: number of measures
- del_m: interval between measures
- lambda: step used in external influence changes

---------Output files----------

---> "ext_corr.dat"
- odd lines contain the correlation between nodes with state 0 and 
their n0 influence
- even lines contain the correlation between nodes with state 1 and 
their n1 influence

---> "network.dat"
- considering the network the vector state which represents the network, 
we count how many sequential nodes have the same state
- odd lines contain number of sequential same state nodes
- even lines contain the correspondend state 

---> "ext_sum.dat"
- first column contains total sum of external influence N0 in the network
- second column contains total sum of external influence N1 in the network
