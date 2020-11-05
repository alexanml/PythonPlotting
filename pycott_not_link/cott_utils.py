"""Module for getting information from COTT input/output files."""
import numpy as np
import knuplot


class TimeZoneException(Exception):
    """Exception for violation of valid time-zones"""
    pass


class MissingInpVarException(Exception):
    """Exception for not finding a variable in a COTT input file"""
    pass


class InvalidInpVarException(Exception):
    """Exception for invalid or unsupported field in COTT input file"""
    pass


def load_tec(filepath):
    """Load a tec-file into a knuplot Tecplot-object.

    Arguments
        filepath    Path to .tec-file
    Returns
        A knuplot Tecplot-object
    """
    return knuplot.Tecplot(filepath)


def get_time(izone, tec):
    """Get the time of a given zone

    Arguments
        izone   Index of time-zone in tec-file.
        tec     Knuplot Tecplot object
    """
    if izone < len(tec.ts):
        return tec.ts[izone]["t"][0]
    else:
        raise TimeZoneException(
            "Tec-object does not have that many time zones.")


def get_last_zone(time, tec):
    """Get the last zone index before a given time.

    If time is negative, return final index.

    Arguments
        time    Time to get last zone index before (s)
        tec     Knuplot Tecplot object
    Returns
        izone   Index of last time zone before given time.
    """
    nzones = len(tec.ts)

    if time >= 0:
        for i in xrange(nzones):
            t = tec.ts[i]["t"][0]
            if t > time:
                izone = i-1
                break
        else:
            raise TimeZoneException(
                "Time in tec-object does not go beyond t=%.3fs" % time)
    else:
        izone = nzones-1
    return izone


def get_closest_zone(time, tec):
    """Get the closest zone index to a given time.

    Arguments
        time    Time (s)
        tec     Knuplot Tecplot object
    Returns
        izone   Index of time zone closest to given time.
    """
    bestdiff = 1e10
    bestzone = -1
    for i in range(0, len(tec.ts)):
        zonetime = get_time(i, tec)
        if abs(zonetime - time) < bestdiff:
            bestdiff = abs(zonetime - time)
            bestzone = i
    return bestzone


def get_var_vs_t(varstr, icell, tec):
    """Get times and a variable at these times, from a given cell.

    Arguments
        varstr      Name of variable in tec-object.
        icell       Index of cell to get data from.
        tec         Knuplot Tecplot object
    Returns
        t           Array of times (s)
        var         Array of variable values (unknown unit)
    """
    t = []
    var = []
    for zonedata in tec.ts:
        tval = zonedata["t"][0]
        varval = zonedata[varstr][icell]
        t.append(tval)
        var.append(varval)
    t = np.array(t)
    var = np.array(var)
    return t, var

def get_var_vs_t_at_x(varstr, x, tec):
    """
    """

    x_from_tec = tec.ts[0]['x']
    icell = (np.abs(x_from_tec - x)).argmin()
    actual_x = x_from_tec[icell]

    t = []
    var = []
    for zonedata in tec.ts:
        tval = zonedata["t"][0]
        varval = zonedata[varstr][icell]
        t.append(tval)
        var.append(varval)
    t = np.array(t)
    var = np.array(var)
    return t, var, actual_x

def get_var_vs_x(varstr, t, tec):
    """Get the x positions and the variable's values at time t
    """

    x = tec.ts[0]['x']
    # get the time zone index closest to the time you want to check
    it = get_closest_zone(t, tec)
    # get the variable at the time with index it for all x-positions
    var = tec.ts[it][varstr]
    x = np.array(x)
    var = np.array(var)
    return x, var

def get_vars_vs_x(varstr_list, t_list, tec):
    """Get the x positions and the variable(s)' values over all the time(s)
    
    Arguments
        varstr_list    Names of variable in tec-object (list of strings)
        t_list         Times at which the variable(s)' values over x should be found (list)
        tec       Knuplot Tecplot object
    Returns
        x         Array of x-position
        var_vals       Array of variable(s)' values over x at the different times

        (The shape of var_vals is n_times x n_vars, where n_times is number of times in
        t_list and n_vars is the number of variables in varstr_list. Each element is
        a list of a variable's values over all x at a specific time.)
    """
    n_vars = len(varstr_list)
    n_times = len(t_list)

    x = tec.ts[0]['x']
    var_vals = np.zeros((n_times, n_vars, len(x)))
    for i in range(n_vars):
        for j in range(n_times):
            # get the time zone index closest to the time you want to check
            it = get_closest_zone(t_list[j], tec)
            # get the variable at the time with index it for all x-positions
            var_vals[j, i, :] = tec.ts[it][varstr_list[i]]
    return x, var_vals

def get_lowest_value(varstr, tec):
    """Get the lowest value occurring of a variable

    Searches for lowest value in both space and time.

    Arguments
        varstr      Name of variable in tec-object.
        tec         Knuplot Tecplot object
    Returns
        varmin      Lowest value found
        icell_min   Cell index where it was found
        izone_min   Time zone where it was found
    """
    varmin = float("inf")

    for izone, zone in enumerate(tec.ts):
        vals = zone[varstr]
        thismin = np.min(vals)
        if thismin < varmin:
            varmin = thismin
            icell_min = np.argmin(vals)
            izone_min = izone
    return varmin, icell_min, izone_min


def get_lowest_value_cell(varstr, icell, tec):
    """Get the lowest value occurring of a variable at given cell

    Searches for lowest value in time.

    Arguments
        varstr      Name of variable in tec-object.
        icell       Index of cell to get data from.
        tec         Knuplot Tecplot object
    Returns
        varmin      Lowest value found
        izone_min   Time zone where it was found
    """
    varmin = float("inf")

    for izone, zone in enumerate(tec.ts):
        vals = zone[varstr]
        thismin = vals[icell]
        if thismin < varmin:
            varmin = thismin
            izone_min = izone
    return varmin, izone_min


def set_userinp_value(filepath, field, value):
    """Set a value in the given COTT input file.

    Arguments
        filepath    Path to user.inp file.
        field       A field name
        value       Value to set
    """
    ui = UserInp(filepath)
    ui.set(field, value)
    ui.write(filepath)


def get_userinp_values(filepath, fields):
    """Get values from a COTT input-file.

    Arguments
        filepath    Path to user.inp file.
        fields      List of field names.
    Returns
        vals        List of strings corresponding to "fields".
    """
    vals = []
    ui = UserInp(filepath)
    for field in fields:
        try:
            vals.append(ui.get(field))
        except AttributeError:
            raise MissingInpVarException("Could not find field %s in %s"
                                         % (field, filepath))
    return tuple(vals)


class UserInp(object):
    """An object to write and read from user.inp files.

    Initialized from a file or as an empty object.
    Set and retrieve values using UserInp.set and UserInp.get.
    To write the object to a user.inp file, call UserInp.write.
    """
    def __init__(self, filename=None):
        """Initializes from given filename.
        If no filename is given, the object is initially empty."""
        if filename:
            # Open file with a+ to allow opening nonexisting files
            with open(filename, "a+") as f:
                f.seek(0)
                self.lines = [line.strip() for line in f.readlines()
                              if len(line.strip()) > 0]
        else:
            self.lines = []

    def set(self, name, value):
        """Set a value.

        This does not write to file, which should be done by calling
        UserInp.write afterwards."""
        lines = self.lines
        linefound = False
        for i, line in enumerate(lines):
            if line.startswith(name + "="):
                line = name + "=" + str(value)
                linefound = True
            if line.startswith("expr." + name + "="):
                line = "expr." + name + "=" + str(value)
            lines[i] = line
        if not linefound:
            lines.append(name + "=" + str(value))

    def get(self, name):
        """Get a value."""
        lines = self.lines
        for line in lines:
            if line.startswith(name + "="):
                return line.split("=")[-1].strip()
        raise AttributeError("Could not find field " + name)

    def write(self, filename):
        """Write to given user.inp file."""
        with open(filename, "w") as f:
            for line in self.lines:
                f.write(line + "\n")

    def category(self, name):
        """Return a UserInpCategory with the given name
        ('spi.num', 'spi.cas' etc.)
        """
        return UserInpCategory(self, name)


class UserInpCategory(object):
    """Represents a category in user.inp.

    This allows setting and getting without specifying spi.(cat). as prefix.
    """
    def __init__(self, userinp, category):
        """Initialize with given UserInp object and a category
        ('spi.cas', 'spi.num') etc.
        """
        self.userinp = userinp
        self.category = category

    def set(self, *args, **kwargs):
        """Set values either using set(name, value)
        or set(name1=value1, name2=value2, ...).
        """
        if len(args) > 0:
            self.userinp.set(self._fullname(args[0]), args[1])
        if len(kwargs) > 0:
            for name, value in kwargs.iteritems():
                self.userinp.set(self._fullname(name), value)

    def get(self, name):
        return self.userinp.get(self._fullname(name))

    def _fullname(self, name):
        """Concatenates name of category and name of field.

        Example:
        UserInpCategory(userinp, 'foo.bar')._fullname('baz') == 'foo.bar.baz'
        """
        return self.category + '.' + name


def get_num_comp(filepath):
    """ Get the number of components from a COTT input-file.

    Arguments
        filepath    Path to user.inp file.
    Returns
        nc          Number of components
    """
    return int(get_userinp_values(filepath, ["spi.cas.nc"])[0])


def get_thermopack_settings(filepath):
    """Get thermopack settings from a COTT input-file.

    Return tuple can be send directly into pytp.tp.init()

    Arguments
        filepath    Path to user.inp file.
    Returns
        eoslib      String for EoS-library
        eos         String for EoS
        mixrule     String for mixing rule
        alpha       String for alpha expression (always Classic)
        nc          Number of components
        compstring  String of comma separated component names
        nphases     Maximum number of phases (always nc+1)

    """
    nc = get_num_comp(filepath)
    #if nc == 1:
    if False:
        ltp, compstring, eos, mixrule = get_userinp_values(
            filepath,
            ("spi.scub.thermopack",
             "spi.scub.comp",
             "spi.scub.ieos",
             "spi.scub.imix")
            )
        ltp = bool(int(ltp))
        if ltp:
            eoslib = "Thermopack"
        else:
            eoslib = "TPlib"
            #raise InvalidInpVarException("COTT not using Thermopack")

    else:
        eoslib, eos, mixrule, compstring = get_userinp_values(
            filepath,
            ("spi.eos.ieoslib",
             "spi.eos.ieos",
             "spi.eos.imix",
             "spi.eos.components")
            )
    alpha = "Classic"
    nphases = 3
    nc_str = len(compstring.strip(" ,").split(","))
    assert nc == nc_str
    return eoslib, eos, mixrule, alpha, nc, compstring, nphases


def get_molar_composition(filepath):
    """Get molar composition from COTT input-file.

    Arguments
        filepath    Path to user.inp file.
    Returns
        z           Array of molar compositions.

    """
    nc = get_num_comp(filepath)
    if nc == 1:
        z = np.array([1.0])
    else:
        fracunit, z_str = get_userinp_values(filepath, ("spi.eos.fracunit",
                                                        "spi.eos.fractions"))
        if fracunit.strip() == "mol":
            z = z_str.strip(" ,").split(",")
            z = [float(zval) for zval in z]
            z = np.array(z)
        else:
            raise NotImplementedError(
                "Mass fractions in inp-file." +
                " Conversion to molar not implemented.")
    return z
