Quick overview of what is included
1. generate_BLM.py --> for generating grid lattice of 4 site lipid model in a bilayer form (1-head, 2-intermediate, 3-tail, 3-tail)
2. config files --> most used for this paper
     - PIEZO CG in vacuum
     - PIEZO CG in mid size bilayer (45k lipids)
     - small size bilayer (20k lipids)
     - mid size bilayer (45k lipids)
     - large size bilayer (180k lipids)
3. Control Files for new parameter (fully annotated) --> perfect for full understanding modification
4. Sample test run on PIEZO CG in vacuum (all files included) tutorial
     - PIEZO_CG_in_vacuum.config
     - production1.control
     - restart_files/ -- (empty)
     - sample_submission_script.sh (specific for stampede3)
5. Kc calculation
     - Spectrum
     - Spectrum run submission sample
     - Kc calculation from output and graphing (ipynb)
6. Tension
     - How to grab data from output
     - graphing (ipynb)
7. PIEZO footprint
     - README included
     - c++ files
8. END
