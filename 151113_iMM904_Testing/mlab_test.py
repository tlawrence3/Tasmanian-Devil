from pymatbridge import Matlab
mlab = Matlab()
mlab.start()
a = 'help_today'
mlab.get_variable('a')

%load_ext pymatbridge
%%matlab
a

mlab.stop()
