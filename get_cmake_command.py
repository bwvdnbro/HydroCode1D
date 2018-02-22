#! /usr/bin/python

################################################################################
# This file is part of HydroCode1D
# Copyright (C) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# HydroCode1D is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HydroCode1D is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with HydroCode1D. If not, see <http://www.gnu.org/licenses/>.
################################################################################

##
# @file get_cmake_command.py
#
# @brief Script that generates the cmake command necessary to configure the code
# with a specific configuration.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

# default options: do not touch this!

configuration_options = {
"rmin_in_m": 0.,
"rmax_in_m": 1.,
"ncell": 1000,
"gamma": 5./3.,
"maxtime_in_s": 1.,
"number_of_snaps": 1,
"eos": "EOS_IDEAL",
"boundaries": "BOUNDARIES_OPEN",
"isothermal_temperature_in_k": 500.,
"potential": "POTENTIAL_POINT_MASS",
"g_internal": 1.,
"mass_point_mass_in_kg": 1.,
"courant_factor": 0.4,
"riemannsolver_type": "RIEMANNSOLVER_TYPE_EXACT",
"dimensionality": "DIMENSIONALITY_1D",
"hydro_order": "HYDRO_ORDER_2",
}

##
# @brief Generate the cmake command to configure the code with a specific
# configuration.
#
# @param custom_options Configuration options that should replace their
# respective default values.
# @param folder Folder where the CMakeLists.txt file is located.
# @return CMake configuration command.
##
def get_cmake_command(custom_options = {}, folder = ".."):
  global configuration_options

  command = "cmake -DCMAKE_BUILD_TYPE=Release"
  for option in configuration_options:
    if option in custom_options:
      value = custom_options[option]
      custom_options[option] = "read"
    else:
      value = configuration_options[option]
    command += " -D{0}={1}".format(option, value)
  command += " " + folder

  # check that all custom options were actually used
  for option in custom_options:
    if not custom_options[option] == "read":
      print "Unknown option:", option
      exit()

  return command

##
# @brief Get the default CMake configuration command.
##
if __name__ == "__main__":
  print get_cmake_command()
