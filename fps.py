import os
import json
from collections import OrderedDict
import numpy as np
import pandas as pd
from . import av_functions


package_directory = os.path.dirname(__file__)
dye_definition = json.load(open(os.path.join(package_directory, 'dye_definition.json')))
dye_names = dye_definition.keys()


class AV(object):
    @property
    def ng(self):
        """The number of grid-points in each direction
        """
        return self.density.shape[0]

    @property
    def density(self):
        """A three-dimensional array of the accessible volume. The values of the array are bigger than one
        if the dye can reach the point.
        """
        return self._density

    @density.setter
    def density(self, v):
        self._density = v

    @property
    def Rmp(self):
        """The mean position of the accessible volume (average x, y, z coordinate)
        """
        return self.points.mean(axis=0)

    def dRmp(self, av):
        """Calculate the distance between the mean positions with respect to the accessible volume `av`

        Remark: The distance between the mean positions is not a real distance

        :param av: accessible volume object
        :return:

        Examples
        --------

        """
        return av_functions.dRmp(self, av)

    def dRDA(self, av):
        """Calculate the mean distance to the second accessible volume (as observed by TCSPC-measurements)

        :param av:
        :return:

        Examples
        --------

        """
        return av_functions.RDAMean(self, av)

    def dRDAE(self, av):
        """Calculate the FRET-averaged mean distance to the second accessible volume (as observed by intensity
        based FRET-measurements)

        :param av:
        :return:

        Examples
        --------

        """
        return av_functions.RDAMeanE(self, av)

    def save_xyz(self, filename):
        """Saves the accessible volume as xyz-file

        :param filename: string
        """
        points = self.points
        fp = open(filename, 'w')
        npoints = len(points)
        fp.write('%i\n' % npoints)
        fp.write('Name\n')
        for p in points:
            fp.write('D %.3f %.3f %.3f\n' % (p[0], p[1], p[2]))
        fp.close()

    def __init__(self, points, density, x0, parameters):
        self.points = points
        self._density = density
        self.x0 = x0
        self.parameters = parameters

    def __len__(self):
        return self.points.shape[0]

    def __str__(self):
        s = 'Acessible Volume\n'
        s += '----------------\n\n'
        s += 'n-points   : %i\n' % len(self)
        s += 'attachment : %.3f, %.3f, %.3f\n\n' % (self.x0[0], self.x0[1], self.x0[2])
        s += 'Parameters: \n'
        for key in self.parameters.keys():
            s += "\t%s: \t %s\n" % (key, self.parameters[key])
        return s


class AVTrajectory(object):
    """
    Examples
    --------

    >>> import mdtraj as md
    >>> import mdtraj.fluorescence

    Load some trajectory

    >>> traj = md.load('hGBP1_out_3.h5')

    Make a new accessible volume trajectory. This places a dye on the specified atom. The dye parameters are either
    passes as or taken from a dye-library (see dye_definition.json). If no dye-parameters are passed default
    parameters are used (not recommended).

    >>> av_traj_1 = mdtraj.fluorescence.fps.AVTrajectory(traj, '18D', attachment_atom_selection='resSeq 7 and name CB')

    For visual inspection the accessible volume can be saved as xyz-file.

    >>> av_traj[0].save_xyz('test_344.xyz')

    A preset dye-parameter is loaded using the argument `dye_parameter_set`. Here the string has to match the string in the
    dye-definition.json file.

    >>> av_traj_2 = mdtraj.fluorescence.fps.AVTrajectory(traj, '18D', attachment_atom_selection='resSeq 7 and name CB', dye_parameter_set='D3Alexa488')

    """

    def __init__(self, trajectory, position_name, dye_parameter_set=None, **kwargs):
        # Trajectory
        self.verbose = kwargs.get('verbose', False)
        self.trajectory = trajectory
        self._avs = dict()

        # Determine vdw-radii from topology
        self.vdw = kwargs.get('vdw', av_functions.get_vdw(trajectory))

        # AV-position
        self.position_name = position_name
        self.simulation_type = kwargs.get('simulation_type', 'AV1')
        attachment_atom_index = kwargs.get('attachment_atom_index', None)
        attachment_atom_selection = kwargs.get('attachment_atom_selection', None)

        if self.verbose:
            print "Attachment atom"
            print attachment_atom_selection
            print type(attachment_atom_selection)

        # Determine attachment atom index
        if isinstance(attachment_atom_index, int):
            self.attachment_atom_index = attachment_atom_index
        elif isinstance(attachment_atom_selection, unicode) or isinstance(attachment_atom_selection, str):
            topology = trajectory.topology
            index = topology.select(attachment_atom_selection)
            if self.verbose:
                print "Using selection string to determine attachment position"
                print "Attachment atom: %s" % index
            if len(index) != 1:
                raise ValueError("The selection does not result in a single atom. Please change your selection.")
            else:
                self.attachment_atom_index = index[0]
        else:
            raise ValueError("Provide either the attachment atom index or a selection string.")

        # Determine either the AV-parameters by the preset dye-parameter
        # or use the parameters provided by the user
        if isinstance(dye_parameter_set, str):
            p = dye_definition[dye_parameter_set]
            self.linker_length = p['linker_length']
            self.linker_width = p['linker_width']
            self.radius_1 = p['radius1']
            self.radius_2 = p['radius2']
            self.radius_3 = p['radius3']
            self.allowed_sphere_radius = p['allowed_sphere_radius']
            self.simulation_grid_resolution = p['simulation_grid_resolution']
        else:
            self.linker_length = kwargs.get('linker_length', 20.0)
            self.linker_width = kwargs.get('linker_width', 0.5)
            self.radius_1 = kwargs.get('radius_1', 4.5)
            self.radius_2 = kwargs.get('radius_2', 4.5)
            self.radius_3 = kwargs.get('radius_3', 4.5)
            self.allowed_sphere_radius = kwargs.get('allowed_sphere_radius', 1.0)
            self.simulation_grid_resolution = kwargs.get('simulation_grid_resolution', 0.5)


    @property
    def _sim_parameter(self):
        re = {
            'vdw': self.vdw,
            'l': self.linker_length,
            'w': self.linker_width,
            'r1': self.radius_1,
            'r2': self.radius_1,
            'r3': self.radius_2,
            'linkersphere': self.allowed_sphere_radius,
            'dg': self.simulation_grid_resolution,
            'atom_i': self.attachment_atom_index
        }
        return re

    def __len__(self):
        return len(self.trajectory)

    def __getitem__(self, key):
        simulation_type = self.simulation_type

        if isinstance(key, int):
            frame_idx = [key]
        else:
            start = 0 if key.start is None else key.start
            stop = None if key.stop is None else key.stop
            step = 1 if key.step is None else key.step
            frame_idx = range(start, min(stop, len(self)), step)

        re = []
        for frame_i in frame_idx:
            # TODO: BETTER conversion nano-meter angstrom
            frame = self.trajectory[frame_i]
            x = (frame.xyz[0, :, 0] * 10.0).astype(np.float64, order='C')
            y = (frame.xyz[0, :, 1] * 10.0).astype(np.float64, order='C')
            z = (frame.xyz[0, :, 2] * 10.0).astype(np.float64, order='C')

            # if av was already calculated use pre-calculated av
            if frame_i in self._avs.keys():
                av = self._avs[frame_i]
            else:
                parameters = self._sim_parameter
                if simulation_type == 'AV1':
                    points, density, ng, x0 = av_functions.calculate_1_radius(x, y, z, **parameters)
                elif simulation_type == 'AV3':
                    points, density, ng, x0 = av_functions.calculate_3_radius(x, y, z, **parameters)
                av = AV(points, density, x0, parameters)
                self._avs[frame_i] = av
            re.append(av)

        if len(re) == 1:
            return re[0]
        else:
            return re


class AvDistanceTrajectory(object):
    """
    The AvPotential class provides the possibility to calculate the reduced or unreduced chi2 given a set of
    labeling positions and experimental distances. Here the labeling positions and distances are provided as
    dictionaries.

    Examples
    --------

    distance_file = './example/hGBP1_distance.json'
    av_dist = mdtraj.fluorescence.fps.AvDistanceTrajectory(traj, distance_file)
    av_dist[:3]
    """

    def __init__(self, trajectory, distance_file, **kwargs):
        d = json.load(open(distance_file, 'r'))
        self.distances = d['Distances']
        self.positions = d['Positions']
        self.n_av_samples = kwargs.get('av_samples', 10000)
        self.bins = kwargs.get('hist_bins', np.linspace(10, 100, 90))
        self.trajectory = trajectory
        self._d = dict()
        self.verbose = kwargs.get('verbose', False)
        self.vdw = av_functions.get_vdw(trajectory)

    def get_avs(self, traj_index):
        frame = self.trajectory[traj_index]
        re = OrderedDict()
        arguments = [
            dict(
                {
                    'vdw': self.vdw,
                    'trajectory': frame,
                    'position_name': position_key,
                },
                **self.positions[position_key]
            )
            for position_key in self.positions
        ]
        avs = map(lambda x: AVTrajectory(**x), arguments)
        for i, position_key in enumerate(self.positions):
            re[position_key] = avs[i]
        return re

    def __len__(self):
        return len(self.trajectory)

    def __getitem__(self, key):
        if isinstance(key, int):
            frame_idx = [key]
        else:
            start = 0 if key.start is None else key.start
            stop = None if key.stop is None else key.stop
            step = 1 if key.step is None else key.step
            frame_idx = range(start, min(stop, len(self)), step)

        re = dict((key, {'rMP': [], 'rDA': [], 'rDAE': [], 'pRDA': [], 'chi2': []}) for key in self.distances.keys())
        for frame_i in frame_idx:
            # Don't repeat calculations
            if frame_i in self._d.keys():
                rDA, rDAE, rMP, chi2, pRDA = self._d[frame_i]
            else:
                # calculate the AVs of the frame
                avs = self.get_avs(frame_i)
                # Calculate the distances
                for distance_key in self.distances:
                    distance = self.distances[distance_key]
                    av1 = avs[distance['position1_name']][0]
                    av2 = avs[distance['position2_name']][0]
                    R0 = distance['Forster_radius']

                    ran_dist = av_functions.random_distances(av1, av2, self.n_av_samples)
                    bin_edges, pRDA = av_functions.pRDA(distances=ran_dist, bins=self.bins, n_samples=self.n_av_samples)
                    rDA = av_functions.RDAMean(distances=ran_dist, nSamples=self.n_av_samples)
                    rDAE = av_functions.RDAMeanE(distances=ran_dist, forster_radius=R0, nSamples=self.n_av_samples)
                    rMP = av_functions.dRmp(av1, av2)
                    if self.verbose:
                        print "RDA: %s" % rDA
                        print "RDA_E: %s" % rDAE
                        print "RDA_mp: %s" % rMP
                    if self.distances[distance_key]['distance_type'] == 'RDAMean':
                        self.distances[distance_key]['model_distance'] = rDA
                    elif self.distances[distance_key]['distance_type'] == 'RDAMeanE':
                        self.distances[distance_key]['model_distance'] = rDAE
                    elif self.distances[distance_key]['distance_type'] == 'Rmp':
                        self.distances[distance_key]['model_distance'] = rMP

                # compare to experiment: calculate chi2
                chi2 = 0.0
                for distance in list(self.distances.values()):
                    dm = distance['model_distance']
                    de = distance['distance']
                    error_neg = distance['error_neg']
                    error_pos = distance['error_pos']
                    d = dm - de
                    chi2 += (d / error_neg) ** 2 if d < 0 else (d / error_pos) ** 2
                self._d[frame_i] = (rDA, rDAE, rMP, chi2, list(pRDA))

            for distance_key in self.distances:
                re[distance_key]['rMP'].append(rMP)
                re[distance_key]['rDA'].append(rDA)
                re[distance_key]['rDAE'].append(rDAE)
                re[distance_key]['pRDA'].append(pRDA)
                re[distance_key]['chi2'].append(chi2)

        return re