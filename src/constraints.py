'''
constraints
'''

class DensityConstraint(object):
    def __init__(self, volume_frac, density_min = 0, density_max = 1.0):
        self.volfrac = volume_frac
        self.xmin = density_min
        self.xmax = density_max

    def volume_frac(self):
        return self.volfrac

    def density_min(self):
        return self.xmin

    def density_max(self):
        return self.xmax

