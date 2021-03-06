# Copyright (C) 2012 Sergio Rossell
#
# This script is part of the EXAMO software
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>
# 
"""
"""


########################################
# UTILITY FUNCTIONS

def importPickle(fileName, mode = 'rb'):
    """
    Imports a pickle object. By default it reads as binary (rb). Setting
    mode allows other ways of opening the file (e.g. mode = 'r')
    """
    import cPickle as pickle
    f = open(fileName, mode)
    return pickle.load(f)

def exportPickle(obj, fileName, mode = 'wb', protocol = -1):
    """
    Exports an object as a pickle file. By default it writes as binary (wb).
    Setting mode allows other ways of opening the file (e.g. mode = 'w')
    """
    import cPickle as pickle
    f = open(fileName, mode)
    pickle.dump(obj, f, protocol = -1)
    f.close()
