#!/usr/bin/python
import knuplot
import sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
from matplotlib.ticker import ScalarFormatter

# Read arguments
try:
  tecfile = sys.argv[1]
  var = sys.argv[2]
  fps = float(sys.argv[3])
  if len(sys.argv)==5:
    savefile = bool(sys.argv[4])
  else:
    savefile = False
except:
  print 'usage: "python tec_to_movie.py [tecfile] "[variable(s)]" [FPS] [savefile = False]"'
  exit()

varlist = var.split()

# Read TEC file
tec = knuplot.Tecplot(tecfile)
n_times = len(tec.ts)
xvals = tec.ts[0]["x"][1:-1]
times = np.array([ tec.ts[zone]["t"][0] for zone in range(n_times)])

var_arrays = []
for var in varlist:
  var_arrays.append([ tec.ts[zone][var][1:-1] for zone in range(n_times)])


if len(varlist)==1:

  # Set plot limits
  minval = 1e20
  maxval = 1e-20
  for arrays in var_arrays:
    for array in arrays:
      local_max = np.max(array)
      local_min = np.min(array)
      if local_max > maxval:
        maxval = local_max
      if local_min < minval:
        minval = local_min

  # Set up plot
  fig = plt.figure()
  ax = plt.axes(xlim=(xvals[0], xvals[-1]), ylim=(minval, maxval))

  lines = []
  line, = ax.plot([], [], lw=1)
  plt.xlabel("x")
  plt.ylabel(varlist[0])
  plt.grid()

  formatter = ScalarFormatter(useOffset=False)
  ax.yaxis.set_major_formatter(formatter)

  def init():
      line.set_data([], [])
      return line,

  def animate(i):
      line.set_data(xvals, var_arrays[0][i])
      if i==0:
        print "i=%5s:  t=%10.4f ms" % (i,times[i]*1000)
      else:
        print "i=%5s:  t=%10.4f ms,  dt=%10.4f ms" % (i,times[i]*1000,(times[i]-times[i-1])*1000)

      return line,

  # Make animation
  anim = animation.FuncAnimation(fig, animate, init_func=init,
                                 frames=n_times, interval=1000.0/fps, blit=True)
  if savefile:
    try:
      print "Making movie file"
      anim.save('animation.mp4', fps=fps)
    except:
      print "Failed"
  print "Showing animation"
  plt.show()

elif len(varlist)==2:
  
  # Set plot limits
  minvals = []
  maxvals= []
  
  for i in range(2): 
    minvals.append(1e20)
    maxvals.append(1e-20)
    for array in var_arrays[i]:
      local_max = np.max(array)
      local_min = np.min(array)
      if local_max > maxvals[i]:
        maxvals[i] = local_max
      if local_min < minvals[i]:
        minvals[i] = local_min
  

  host = host_subplot(111, axes_class=AA.Axes)
  plt.subplots_adjust(right=0.75)
  par1 = host.twinx()


  host.set_xlim(xvals[0], xvals[-1])
  host.set_ylim(minvals[0], maxvals[0])

  host.set_xlabel("Dist")
  host.set_ylabel(varlist[0])
  par1.set_ylabel(varlist[1])
  #par2.set_ylabel("Velocity")

  p1, = host.plot(xvals, var_arrays[0][-1])
  p2, = par1.plot(xvals, var_arrays[1][-1])

  par1.set_ylim(minvals[1], maxvals[1])

  host.axis["left"].label.set_color(p1.get_color())
  par1.axis["right"].label.set_color(p2.get_color())
  


  plt.show()


