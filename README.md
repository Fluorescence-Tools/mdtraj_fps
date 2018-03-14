# MDTrajFPS
MDTrajFPS is tool to calculate FRET observables from MD-trajectories.

MDTrajFPS depends on MDTraj and utilizes the original libraries published together with the "FRET Positioning System" (FPS) in Nature Methods 9, 1218–1225 (2012) "A toolkit and benchmark study for FRET-restrained high-precision structural modeling". 

The implicit model of the label is being improved to support the latest developments 
[![DOI for Citing COSB](https://img.shields.io/badge/DOI-10.1016/j.sbi.2016.11.012-blue.svg)](https://doi.org/10.1016/j.sbi.2016.11.012). Hence, the library for the calculation of the implicit labels are currently only provided as pre-compiled binaries for (Windows x32, x64, MacOS, and Linux x32, x64). 

With MDTrajFPS, you can

- Process every MD format supported by MDTraj
- Calculate FRET efficienes using a fast implicit dye model
- Save calculated positional distributions of dyes for later vizualization

The original FPS-software package is available at http://www.mpc.hhu.de/software/fps.html.


##  Code Example

```python
import mdtraj as md
import fps

traj = md.load('./exsample/hGBP1_out_3.h5')
av_traj = fps.AVTrajectory(traj, '18D', attachment_atom_selection='resSeq 7 and name CB')

# save accessible volume as xyz-file
av_traj[0].save_xyz('test_344.xyz')

# Use certain dye-parameter set (see dye_defenitiion.json file)
av_traj = fps.AVTrajectory(traj, '18D', attachment_atom_selection='resSeq 7 and name CB', dye_parameter_set='D3Alexa488')

# Calculate
distance_file = './exsample/hGBP1_distance.json'
av_dist = fps.AvDistanceTrajectory(traj, distance_file)
av_dist[0]
```

## Dependencies
The windows dll requires the following redistributables (not included)

MS-visual studio 2005
	vcredist_x64
	vcredist_x86


##  Citations 
* MDTraj - [![DOI for Citing MDTraj](https://img.shields.io/badge/DOI-10.1016%2Fj.bpj.2015.08.015-blue.svg)](http://doi.org/10.1016/j.bpj.2015.08.015)
* FPS - [![DOI for Citing FPS](https://img.shields.io/badge/DOI-10.1038/nmeth.2222-blue.svg)](http://doi.org/10.1038/nmeth.2222)


## License

GNU LGPL version 2.1, or at your option a later version of the license.
Various sub-portions of this library may be independently distributed under
different licenses. See those files for their specific terms.
