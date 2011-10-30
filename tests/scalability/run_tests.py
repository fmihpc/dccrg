#! /usr/bin/env python
# -*- coding: utf-8 -*-

'''
A program for starting general scalability tests and gathering results.

Copyright 2011 Finnish Meteorological Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License version 3
as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
'''

from optparse import OptionParser
from os import chdir, getcwd
from subprocess import check_call


def test_using_mpirun(processes, data_size, solution_time, test_dir_name, options):
	pwd = getcwd()
	chdir(test_dir_name)
	run_command = "mpirun -np " + str(processes)
	run_command += " scalability --data_size " + str(data_size)
	run_command += " --solution_time " + str(solution_time)
	run_command += " --x_length " + str(options.x_length)
	run_command += " --y_length " + str(options.y_length)
	run_command += " --z_length " + str(options.z_length)
	run_command += " --neighborhood_size " + str(options.neighborhood_size)
	run_command += " > logfile"
	check_call(run_command, shell = True)
	chdir(pwd)


parser = OptionParser()

parser.add_option(
	"--prefix",
	type = "string",
	default = "results",
	dest = "prefix",
	help = "Main directory for running the tests",
	metavar = "PREFIX")

parser.add_option(
	"--gather-results",
	default = False,
	action = "store_true",
	dest = "gather_results",
	help = "Gather results of finished tests")

parser.add_option(
	"--process-min",
	type = "int",
	default = 1,
	dest = "process_min",
	help = "Start launching at MINP processes",
	metavar = "MINP")

parser.add_option(
	"--process-max",
	type = "int",
	default = 6,
	dest = "process_max",
	help = "Stop launching at MAXP processes",
	metavar = "MAXP")

parser.add_option(
	"--process-stride",
	type = "int",
	default = 1,
	dest = "process_stride",
	help = "Increment number of launched processes by PSTRIDE",
	metavar = "PSTRIDE")

parser.add_option(
	"--data-size-min",
	type = "int",
	default = 1,
	dest = "data_size_min",
	help = "Start testing at DMIN bytes of data per cell",
	metavar = "DMIN")

parser.add_option(
	"--data-size-max",
	type = "int",
	default = 1000000,
	dest = "data_size_max",
	help = "Stop testing at DMAX bytes of data per cell",
	metavar = "DMAX")

parser.add_option(
	"--data-size-factor",
	type = "float",
	default = 10.0,
	dest = "data_size_factor",
	help = "Multiply data size by DFAC",
	metavar = "DFAC")

parser.add_option(
	"--solution-time-min",
	type = "float",
	default = 1e-6,
	dest = "solution_time_min",
	help = "Start testing at TMIN second solution time per cell",
	metavar = "TMIN")

parser.add_option(
	"--solution-time-max",
	type = "float",
	default = 1.0,
	dest = "solution_time_max",
	help = "Stop testing at TMAX second solution time per cell",
	metavar = "TMAX")

parser.add_option(
	"--solution-time-factor",
	type = "float",
	default = 10.0,
	dest = "solution_time_factor",
	help = "Multiply data solution time by TFAC",
	metavar = "TFAC")

parser.add_option(
	"--x-length",
	type = "int",
	default = 10,
	dest = "x_length",
	help = "Test using a grid with XLEN unrefined cell in x direction",
	metavar = "XLEN")

parser.add_option(
	"--y-length",
	type = "int",
	default = 10,
	dest = "y_length",
	help = "Test using a grid with YLEN unrefined cell in y direction",
	metavar = "YLEN")

parser.add_option(
	"--z-length",
	type = "int",
	default = 10,
	dest = "z_length",
	help = "Test using a grid with ZLEN unrefined cell in z direction",
	metavar = "ZLEN")

parser.add_option(
	"--timesteps",
	type = "int",
	default = 10,
	dest = "timesteps",
	help = "Test with STEPS number of timesteps",
	metavar = "STEPS")

parser.add_option(
	"--neighborhood-size",
	type = "int",
	default = 1,
	dest = "neighborhood_size",
	help = "Test with NHOOD size neighborhood",
	metavar = "NHOOD")

parser.add_option(
	"--launcher-type",
	type = "string",
	default = "mpirun",
	dest = "launcher_type",
	help = "Start test programs with launcher type LTYPE",
	metavar = "LTYPE")

(options, args) = parser.parse_args()

check_call(["mkdir", "-pv", options.prefix])
check_call(["cp", "-v", "scalability", options.prefix])

# create run directories, job scripts, etc and start
for processes in range(options.process_min, options.process_max + 1, options.process_stride):

	data_size = options.data_size_min
	while data_size <= options.data_size_max:

		solution_time = options.solution_time_min
		while solution_time <= options.solution_time_max:

			test_dir_name = options.prefix + "/" + str(processes) + "p_" + str(data_size) + "B_" + str(solution_time) + "s"
			check_call(["mkdir", "-pv", test_dir_name])
			check_call(["ln", "-svf", "../scalability", test_dir_name])

			if options.launcher_type == 'mpirun':
				test_using_mpirun(processes, data_size, solution_time, test_dir_name, options)

			solution_time *= options.solution_time_factor
		data_size = int(data_size * options.data_size_factor)

