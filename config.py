#import pbs.classes.Dynamic
import os
import pbs

#self.execfile(os.path.join('modules/manim/config.py'))
#self.execfile(os.path.join('modules/gplot/config.py'))
#self.execfile(os.path.join('modules/graph/config.py'))
#self.execfile(os.path.join('modules/esolv/config.py'))

l = pbs.LibraryPython(self, 'nbody', __file__)

#l = pbs.classes.Static.Static("graph", self)

#l.make()

#l.l_defines.append('GR_GRAPH_HPP_IN_LOGGER_MODE=logs::mode::RUN_TIME')
#l.l_defines.append('GR_GRAPH_HPP_IN_LOGGER_LEVEL=1')

#l.doc_out_dir = "/media/sf_P_DRIVE/html/graph"

#l.add_dep('esolv')
#l.add_dep('graph')
#l.add_dep('gplot')
#l.add_dep('manim')

#l.args.libraries.append('glpk')
#l.args.libraries.append('python2.7')
#l.args.libraries.append('boost_python-py36')
#l.args.libraries.append('python3.6m')

self.parts.append(l)

#self.execfile(os.path.join(__dir__, "tests/config.py"))

