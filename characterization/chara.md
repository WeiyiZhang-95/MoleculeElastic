# Characterization 

  Quantifying molecules behavior in different condition.

## Contents
  
  - FDcurve:
  
    Featuring molecule behavior while being deformed by changing distance between two ends (in the code: carboxylates in particularly). The calculation relies on mostly forces changing tendency along deformation.
    
  - config folder:
  
    LAMMPS script for optimizing molecules, deforming simulation and hessian data generating.
    
    _\*Noticing that freqSpectrum.py and dump.py only work with python2_
  
    - example:
      _\* The files used are labled as linker{number}+something alse_
      _\* If changes wanted, revise all the codes_
      ```bash
      #!/bin/bash
      DIRECTORY=linker${args[0]}_deformation
      if [ ! -d "$DIRECTORY" ]; then
        # Control will enter here if $DIRECTORY doesn't exist.
              mkdir ${DIRECTORY}
      fi

      cp ${config_dir}/config/* ${DIRECTORY}
      cp linker${args[0]}.coeff ./${DIRECTORY}/
      cp linker${args[0]}.lmpdat ./${DIRECTORY}/
      cd ./${DIRECTORY}/

      ./extract-element-sequence.sh ${args[0]}

      IN=in.deform
      lmp_mpi < $IN > $IN\.out -var LinkerID ${args[0]}
      ## or mpirun lmp_mpi < $IN > $IN\.out -var LinkerID ${args[0]}
      ## lmp_mpi is the excutable file of LAMMPS
      ```

