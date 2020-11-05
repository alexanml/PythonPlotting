# TODO - gui?

"""
Knuplot is a command-driven plotting tool for the tec-producing codes 
like e.g. levelZ and COTT.
"""

# Import necessary modules
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

# This dictionary can be used for labels etc.
# If not present, the datamap will be used instead.
miscdict={'`DA^-' :'Relative area',
          'A^-'   :'Area',
          'u_c'   :'Horizontal velocity',
          'v_c'   :'Vertical velocity',
          'Re_d'  :'Reynolds number',
          't'     :'Time',
          'F_x'   :'Horizontal force',
          'F_y'   :'Vertical force',
          'surfm' :'Surfactant mass'
          }

class Miscplot:
    """Read a misc file"""

    def __init__(self,file='levelZ-misc.tec'):
        self.datamap={}
        self.data=[]
        print("Loading %s..." % file)
        f = open(file)
        # Heading
        f.readline()
        # The time is always outputted
        f.readline()
        idx = 0
        self.datamap[idx] = 't'
        idx = idx + 1
        self.data.append([])
        # Register other variables in datamap
        while True:
            line = f.readline()
            if "ZONE" in line:
                break
            else:
                var = line.split()[-1]
                self.datamap[idx] = eval(var)
                idx = idx + 1
                self.data.append([])
        # Load data (There's probably a "nicer" way of doing this...)
        nvars = idx
        EOF = False
        while not EOF:
            for i in range(0,nvars):
                line = f.readline()
                if not line:
                    EOF = True
                else:
                    try:
                        self.data[i].append(float(line))
                    except:
                        print("Warning, corrupted data for variable %i" % i)
                        print(line)
                        print("Converting to zero.")
                        self.data[i].append(0.0)

        print("%i variables at %i data points loaded!" % (len(self.datamap),
                                                          len(self.data[0])))
        print("The available variables are:")
        for i in range(0,nvars):
            print('%i: %s' % (i,self.datamap[i]))

    def plot(self,y,x=0,fig=-1,style='',shownow=True):
        """Plot a variable in the data array.

        Input arguments:
        y:     The integer index of the variable to plot.
        x:     The x variable. Default is 0, which is time.
        fig:   The figure number to plot on. Default is a new figure.
        style: Plotting style. Default is let pylab choose.
        """
        if fig == -1:
            plt.figure()
        else:
            plt.figure(fig)
            plt.hold(True)
        plot(self.data[x],self.data[y],style)
        try:
            plt.xlabel(miscdict[self.datamap[x]])
        except:
            plt.xlabel(self.datamap[x])
        try:
            plt.ylabel(miscdict[self.datamap[y]])
        except:
            plt.ylabel(self.datamap[y])
        if shownow:
            plt.show()

class Tecplot:
    """I don't know how to write parsers, so this is just a hackish way to read
    tec files"""

    def __init__(self,file='levelZ.tec'):
        print('Loading ' + file)
        f = open(file)
        # Read TITLE
        line = f.readline()
        titleindex = line.find('"')
        endindex = line.find('"', titleindex+1)
        self.title = line[titleindex+1:endindex].strip()
        # Read names of variables from VARIABLES to ZONE
        line = f.readline()
        varnames = line.split('=')[-1]
        line = f.readline()
        while "ZONE" not in line:
            varnames += line
            line = f.readline()
        varnames = varnames.split(',')
        idx = 0
        self.vars = []
        for varname in varnames:
            # Find variable name without "" and whitespace
            varname = varname.replace('"','').split()[0]
            self.vars.append(varname)

        nvars = len(self.vars)
        eof = False
        self.ts = []
        self.zonetitles = []
        zi = 0 # Zone index
        while not eof:
            # Read zone definition
            titleindex = line.find('T=')
            if titleindex > 0:
                endindex = line.find('"', titleindex+3)
                title = line[titleindex + 3:endindex].strip()
                self.zonetitles.append(title)
            else:
                self.zonetitles.append('')
            # Limited support for VARSHARELIST.
            # Support following form: VARSHARELIST=([1-2,5,8])
            # Assume variables should be taken from zone 1
            varsharelistindex = line.find('VARSHARELIST=')
            varshareidxlist = []
            if varsharelistindex > 0:
                startindex = line.find('[', varsharelistindex+1)
                endindex = line.find(']', varsharelistindex+1)
                varsharelist = line[startindex+1:endindex].strip()
                varshares = varsharelist.split(',')
                for varshare in varshares:
                    rangeindex = varshare.find('-')
                    if rangeindex > 0:
                        i0 = int(varshare[0:rangeindex].strip())
                        i1 = int(varshare[rangeindex+1:].strip())
                        for i in range(i1-i0+1):
                            varshareidxlist.append(i0+i-1)
                    else:
                        varshareidxlist.append(int(varshare)-1)

            iindex = line.find('I=')
            jindex = line.find('J=')
            kindex = line.find('K=')
            dim = np.ones(0, dtype='int')
            ndim = 0
            for i in [iindex, jindex, kindex]:
                if i > 0:
                    commaindex = line.find(',', i)
                    dimlength = int(line[i + 2:commaindex])
                    if dimlength > 1:
                        dim = np.append(dim, dimlength)
                        ndim += 1
            self.ts.append({}) # Make a dictionary
            npoint = np.prod(dim)
            for var in self.vars:
                self.ts[zi][var]= np.zeros(npoint)
            for i in range(0, npoint):
                for j, var in enumerate(self.vars):
                    if j in varshareidxlist:
                        self.ts[zi][var][i] = self.ts[0][var][i]
                    else:
                        self.ts[zi][var][i] = float(f.readline())
            for var in self.vars:
                self.ts[zi][var] = np.reshape(self.ts[zi][var], dim, order='F')
            zi += 1
            # Read heading for new zone
            line = f.readline()
            if not line:
                eof = True

        # Number of zones
        self.nzones = zi

        # Print info
        print('Finished loading.')
        print('Grid size:', dim)
        print('Number of zones:',self.nzones)
        print('The available variables are:')
        print(self.vars)

    def zonetotxt(self, filename, zoneindex=0, withBoundary=True):
        """ Writes a single zone to a text file, with space-separated
        columns.
        """
        zone = self.ts[zoneindex]
        if withBoundary:
            zonelen = len(zone[self.vars[0]])
        else:
            zonelen = len(zone[self.vars[0]]) - 2
        matrix = np.zeros((zonelen, len(self.vars)))
        header = '# '
        endIndex = len(self.vars)
        for i in range(0, endIndex):
            var = self.vars[i]
            header += var + ' (' + str(i+1) + ')'
            if i < endIndex-1:
                header += ', '
            if withBoundary:
                matrix[:, i] = zone[var].ravel()
            else:
                matrix[:, i] = zone[var][1:zonelen+1]
        with file(filename, 'w') as outfile:
            outfile.write(header + '\n')
            savetxt(outfile, matrix)

    def plot(self, y, x='x', zone=0, fig=-1, style=''):
        """One-dimensional plot of a variable in the data array.

        Input arguments:
        y:     Name of the variable to plot.
        x:     The x variable. Default is 'x'.
        zone:  Zone index. Default is 0.
        fig:   The figure number to plot on. Default is a new figure.
        style: Plotting style. Default is let pylab choose.
        """
        if fig == -1:
            plt.figure()
        else:
            plt.figure(fig)
            plt.hold(True)
        plt.plot(self.ts[zone][x],self.ts[zone][y],style)
        try:
            plt.xlabel(miscdict[x])
        except:
            plt.xlabel(x)
        try:
            plt.ylabel(miscdict[y])
        except:
            plt.ylabel(y)
        plt.show()

    def surfplot(self,t,var,var2=-1,style='k-',fig=1):
        """Plot a quantity along the zero level set using bilinear
        interpolation.  If two variables are present, a vector is assumed and
        the vector magnitude is plotted"""
        f2=plt.figure(1322)
        d = self.ts[t]
        v = d[var]
        if var2 != -1:
          v2 = d[var2]
        cs = plt.contour(d['x'],d['y'],d['`f'],[0.0])
        plt.close(f2)
        lc = cs.collections[0]
        try:
            # For older versions of matplotlib
            ps = lc.get_verts()
        except:
            # For newer versions
            ps = lc.get_paths()[0].vertices
        npo = len(ps)
        intf=np.zeros(npo)
        arclen=np.zeros(npo)
        k=0
        h = d['x'][2,2] - d['x'][1,2]
        for p in ps:
            xp,yp = p
            i = np.floor(xp/h)
            j = np.floor(yp/h)
            if var2 == -1:
                x1 = (i)*h; x2 = (i+1.)*h
                y1 = (j)*h; y2 = (j+1.)*h
                c11 = v[i,j];   c12 = v[i,j+1]
                c21 = v[i+1,j]; c22 = v[i+1,j+1]
            else:
                x1 = (i)*h; x2 = (i+1)*h
                y1 = (j)*h; y2 = (j+1)*h
                c11 = np.sqrt(v[i,j]**2+v2[i,j]**2)
                c12 = np.sqrt(v[i,j+1]**2+v2[i,j+1]**2)
                c21 = np.sqrt(v[i+1,j]**2+v2[i+1,j]**2)
                c22 = np.sqrt(v[i+1,j+1]**2+v2[i+1,j+1]**2)
            denom = (x2-x1)*(y2-y1)
            intf[k] = (c11*(x2-xp)*(y2-yp)
                       +  c21*(xp-x1)*(y2-yp)
                       +  c12*(x2-xp)*(yp-y1)
                       +  c22*(xp-x1)*(yp-y1))/denom

            if k==0:
                arclen[k]=0;
            else:
                xm,ym = ps[k-1]
                arclen[k]=arclen[k-1] + np.sqrt((xp-xm)**2 + (yp-ym)**2)

            k=k+1

        plt.figure(fig)
        plt.hold(True)
        plt.plot(arclen,intf,style)
        plt.show()
        return arclen,intf

    def addvar(self,varname,func):
        """Add a new variable to the ts array.

        Give as input the name of the new variable, and a function to compute
        it.
        """
        self.vars.append(varname)
        for zone in range(0,self.nzones):
            self.ts[zone][varname]= func(self,zone)

    def slice(self,t,var,pos,dir,style='k-',fig=1):
        """Create a slice plot."""
        d = self.ts[t]
        if dir=='x':
            x = d['x'][:,pos]
            y = d[var][:,pos]
        elif dir=='y':
            x = d['y'][pos,:]
            y = d[var][pos,:]
        else:
            print("No such dir!\n")
        plt.figure(fig)
        plt.plot(x,y,style)
        return x,y

    def cont(self,t,var,vals=[],n=10,fig=1):
        """Create a contour plot."""
        d = self.ts[t]
        if len(vals)>0:
            plt.contour(d['x'],d['y'],d[var],vals)
        else:
            plt.contour(d['x'],d['y'],d[var],n)
        plt.figure(fig)
        plt.hold(True)

    def contf(self,t,var,rng=10,fig=1):
        """Create a contour plot."""
        d = self.ts[t]
        plt.contourf(d['x'],d['y'],d[var],rng,antialiased=True)
        plt.figure(fig)
        plt.hold(True)

    def quiv(self,t,u='u',v='v',fig=1):
        """Create an arrow plot.

        Defaults to velocity plot.
        """
        d = self.ts[t]
        plt.figure(fig)
        plt.quiver(d['x'],d['y'],d[u],d[v])
        plt.show()

    def tecsport(self,filename,varnames,vardata):
        """Export to a Tecplot-loadable file.

        Currently only works for data gotten from slice and surfplot.
        """
        f = open(filename,'w')
        f.write("TITLE=knuplot\n")
        f.write("VARIABLES =\n")
        for var in varnames:
            f.write("\""+var+"\",\n")
        f.write("ZONE T=Globalvar, F=POINT\n")
        n = len(vardata[0])
        for i in range(0,n):
            for j,var in enumerate(varnames):
                f.write(str(vardata[j][i])+"\n")
        f.close()

# Behaviour when run as a script
if __name__=='__main__':
    print('You have to start knuplot from the python shell.')
    print('Type "from knuplot import *", and then start loading data!')
    print('Use Tecplot for levelZ.tec and Miscplot for levelZ-misc.tec')
