'''
constraints
'''


class DensityConstraint(object):
    def __init__(self, volume_frac, density_min=0.0, density_max=1.0, Emin=1e-9):
        self.volfrac = volume_frac
        self.xmin = density_min
        self.xmax = density_max
        self.Emin = Emin

    def volume_frac(self):
        return self.volfrac

    def density_min(self):
        return self.xmin

    def density_max(self):
        return self.xmax

    def zero_stiffness(self):
        return self.Emin

