import numpy as np
from refnx.analysis import Parameter, possibly_create_parameter, Parameters
from refnx.reflect import SLD, Component, Slab

class VolMono(Component):
    def __init__(self, b_heads, thick_heads, b_tails, tanford_length, chain_tilt, molvols,
                 rough_head_tail, rough_preceding_mono, reverse_monolayer=False, name=''):
        super(VolMono, self).__init__()
        self.head_mol_vol = possibly_create_parameter(molvols[0], '{} - head_molecular_volume'.format(name))
        self.tail_mol_vol = possibly_create_parameter(molvols[1], '{} - tail_molecular_volume'.format(name))
        #self.apm = possibly_create_parameter(apm, '{} - area_per_molecule'.format(name))
        if isinstance(b_heads, complex):
            self.b_heads_real = possibly_create_parameter(b_heads.real,
                                                               name='{} - b_heads_real'.format(name))
            self.b_heads_imag = possibly_create_parameter(b_heads.imag,
                                                               name='{} - b_heads_imag'.format(name))
        else:
            self.b_heads_real = possibly_create_parameter(b_heads,
                                                               name='{} - b_heads_real'.format(name))
            self.b_heads_imag = possibly_create_parameter(0, name='{} - b_heads_imag'.format(name))
        if isinstance(b_tails, complex):
            self.b_tails_real = possibly_create_parameter(b_tails.real,
                                                               name='{} - b_tails_real'.format(name))
            self.b_tails_imag = possibly_create_parameter(b_tails.imag,
                                                               name='{} - b_tails_imag'.format(name))
        else:
            self.b_tails_real = possibly_create_parameter(b_tails, name='{} - b_tails_real'.format(name))
            self.b_tails_imag = possibly_create_parameter(0, name='{} - b_tails_imag'.format(name))

        self.thick_heads = possibly_create_parameter(thick_heads, name='{} - thick_heads'.format(name))
        self.tail_length = possibly_create_parameter(tanford_length, name='{} - tail_length'.format(name))
        self.cos_rad_chain_tilt = possibly_create_parameter(chain_tilt, name='{} - chain_tilt'.format(name))

        self.phit = possibly_create_parameter(0., name='{} - phit'.format(name))
        self.phih = possibly_create_parameter(0.5, name='{} - phih'.format(name))

        self.rough_head_tail = possibly_create_parameter(rough_head_tail, name='{} - rough_head_tail'.format(name))
        self.rough_preceding_mono = possibly_create_parameter(rough_preceding_mono,
                                                              name='{} - rough_preceding_mono'.format(name))

        self.solventsld = possibly_create_parameter(10.8, name='{} - solvent sld'.format(name))
        self.solventsldi = possibly_create_parameter(0., name='{} - solvent sldi'.format(name))
        self.solventrough = possibly_create_parameter(3.3, name='{} - solvent rough'.format(name))
        self.supersld = possibly_create_parameter(0., name='{} - super sld'.format(name))
        self.supersldi = possibly_create_parameter(0., name='{} - super sldi'.format(name))

        self.reverse_monolayer = reverse_monolayer
        self.name = name

    @property
    def slabs(self):
        """
        Returns
        -------
        slab_model = array of np.ndarray
            Slab representaions of monolayer
        """
        layers = np.zeros((4, 5))

        layers[0, 0] = 10.
        layers[0, 1] = self.supersld
        layers[0, 2] = self.supersldi
        layers[0, 3] = 0.
        layers[0, 4] = 0.

        layers[1, 0] = self.tail_length * self.cos_rad_chain_tilt
        layers[1, 1] = self.b_tails_real * 1.e6 / self.tail_mol_vol * (1 - self.phit) + self.supersld * (self.phit)
        layers[1, 2] = self.b_tails_imag * 1.e6 / self.tail_mol_vol * (1 - self.phit) + self.supersldi * (self.phit)
        layers[1, 3] = self.rough_preceding_mono
        layers[1, 4] = 0.

        layers[2, 0] = self.thick_heads
        layers[2, 1] = self.b_heads_real * 1.e6 / self.head_mol_vol * (1 - self.phih) + self.solventsld * (self.phih)
        layers[2, 2] = self.b_heads_imag * 1.e6 / self.head_mol_vol * (1 - self.phih) + self.solventsldi * (self.phih)
        layers[2, 3] = self.rough_head_tail
        layers[2, 4] = 0.

        layers[3, 0] = self.solventrough + self.solventrough * 1.5
        layers[3, 1] = self.solventsld
        layers[3, 2] = self.solventsldi
        layers[3, 3] = self.solventrough

        return layers

    @property
    def parameters(self):
        p = Parameters(name=self.name)
        p.extend([self.b_heads_real, self.b_heads_imag, self.b_tails_real, self.b_tails_imag, self.thick_heads,
                  self.tail_length, self.cos_rad_chain_tilt, self.tail_mol_vol, self.head_mol_vol, self.rough_head_tail,
                  self.rough_preceding_mono, self.solventrough, self.phit, self.phih])
        return p

    def lnprob(self):
        return 0


def set_contraints(lipids, vary_tails=False):
    for i in range(1, len(lipids)):
        #lipids[i].thick_heads.constraint = lipids[0].thick_heads
        if not vary_tails:
            lipids[i].tail_length.constraint = lipids[0].tail_length
        #thick_tail = lipids[0].tail_length.value ** 2 * np.cos(np.deg2rad(lipids[0].chain_tilt.value)) ** 2
        #num = thick_tail * lipids[0].apm.value **2 * lipids[i].tail_mol_vol.value
        #den = lipids[i].tail_length.value ** 2 * lipids[i].apm.value ** 2 * lipids[0].tail_mol_vol.value
        #lipids[i].chain_tilt.constraint = np.rad2deg(np.arccos(np.sqrt(num / den)))
        #lipids[i].rough_head_tail.constraint = lipids[0].rough_head_tail
        #lipids[i].rough_preceding_mono.constraint = lipids[0].rough_preceding_mono
        lipids[i].tail_mol_vol.constraint = lipids[0].tail_mol_vol
        lipids[i].head_mol_vol.constraint = lipids[0].head_mol_vol
        #lipids[i].solventrough.constraint = lipids[0].solventrough.constraint
    return lipids
