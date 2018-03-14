import mdtraj as md
import fps

traj = md.load('./examples/hGBP1_out_3.h5')
av_traj = fps.AVTrajectory(traj, '18D', attachment_atom_selection='resSeq 7 and name CB')

# save accessible volume as xyz-file
av_traj[0].save_xyz('test_344.xyz')

# Use certain dye-parameter set (see dye_defenitiion.json file)
av_traj = fps.AVTrajectory(traj, '18D', attachment_atom_selection='resSeq 7 and name CB', dye_parameter_set='D3Alexa488')

# Calculate
distance_file = './examples/hGBP1_distance.json'
av_dist = fps.AvDistanceTrajectory(traj, distance_file)
print av_dist[0]
