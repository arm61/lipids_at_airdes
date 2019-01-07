import numpy as np
import periodictable as pt
from refnx.analysis import possibly_create_parameter, Parameters
from refnx.reflect import Component


class VolMono(Component):
    def __init__(self, b_heads, thick_heads, b_tails, tanford_length, molvols,
        reverse_monolayer=True, name=''):
        super(VolMono, self).__init__()
        self.head_mol_vol = possibly_create_parameter(
                molvols[0], '{} - head_molecular_volume'.format(name))
        self.tail_mol_vol = possibly_create_parameter(
                molvols[1], '{} - tail_molecular_volume'.format(name))
        if isinstance(b_heads, complex):
            self.b_heads_real = possibly_create_parameter(
                    b_heads.real, name='{} - b_heads_real'.format(name))
            self.b_heads_imag = possibly_create_parameter(
                    b_heads.imag, name='{} - b_heads_imag'.format(name))
        else:
            self.b_heads_real = possibly_create_parameter(
                    b_heads, name='{} - b_heads_real'.format(name))
            self.b_heads_imag = possibly_create_parameter(
                    0, name='{} - b_heads_imag'.format(name))
        if isinstance(b_tails, complex):
            self.b_tails_real = possibly_create_parameter(
                    b_tails.real, name='{} - b_tails_real'.format(name))
            self.b_tails_imag = possibly_create_parameter(
                    b_tails.imag, name='{} - b_tails_imag'.format(name))
        else:
            self.b_tails_real = possibly_create_parameter(
                    b_tails, name='{} - b_tails_real'.format(name))
            self.b_tails_imag = possibly_create_parameter(
                    0, name='{} - b_tails_imag'.format(name))

        self.thick_heads = possibly_create_parameter(
                thick_heads, name='{} - thick_heads'.format(name))
        self.thick_tails = possibly_create_parameter(
                tanford_length, name='{} - thick_tails'.format(name))

        self.phit = possibly_create_parameter(
                0., name='{} - phit'.format(name))
        self.phih = possibly_create_parameter(
                0.5, name='{} - phih'.format(name))

        self.rough_head_tail = possibly_create_parameter(
                3., name='{} - rough_head_tail'.format(name))
        self.rough_preceding_mono = possibly_create_parameter(
                3., name='{} - rough_preceding_mono'.format(name))

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
        layers = np.zeros((2, 5))

        layers[0, 0] = self.thick_tails
        layers[0, 1] = self.b_tails_real * 1.e6 / self.tail_mol_vol
        layers[0, 2] = self.b_tails_imag * 1.e6 / self.tail_mol_vol
        layers[0, 3] = self.rough_preceding_mono
        layers[0, 4] = self.phit

        layers[1, 0] = self.thick_heads
        layers[1, 1] = self.b_heads_real * 1.e6 / self.head_mol_vol
        layers[1, 2] = self.b_heads_imag * 1.e6 / self.head_mol_vol
        layers[1, 3] = self.rough_head_tail
        layers[1, 4] = self.phih

        return layers

    @property
    def parameters(self):
        p = Parameters(name=self.name)
        p.extend([self.b_heads_real, self.b_heads_imag, self.b_tails_real,
                  self.b_tails_imag, self.thick_heads, self.thick_tails,
                  self.tail_mol_vol, self.head_mol_vol, self.rough_head_tail,
                  self.rough_preceding_mono, self.phit, self.phih])
        return p

    def lnprob(self):
        return 0


def set_constraints(lipids, vary_tails=False):
    for i in range(1, len(lipids)):
        lipids[i].thick_heads.constraint = lipids[0].thick_heads
        if not vary_tails:
            lipids[i].thick_tails.constraint = lipids[0].thick_tails
        lipids[i].phit.constraint = lipids[0].phit
        lipids[i].tail_mol_vol.constraint = lipids[0].tail_mol_vol
        lipids[i].head_mol_vol.constraint = lipids[0].head_mol_vol
    return lipids


def get_scattering_length(component, neutron=False):
    if not neutron:
        scattering_length = 0 + 0j
        import scipy.constants as const
        cre = const.physical_constants['classical electron radius'][0] * 1e10
        for key in component:
            scattering_length += np.multiply(
                    pt.elements.symbol(key).xray.scattering_factors(
                            energy=12)[0], cre * component[key])
            scattering_length += np.multiply(
                    pt.elements.symbol(key).xray.scattering_factors(
                            energy=12)[1], cre * component[key]) * 1j
    else:
        scattering_length = 0 + 0j
        for key in component:
            scattering_length += (pt.elements.symbol(key).neutron.b_c *
                                  component[key])
            if pt.elements.symbol(key).neutron.b_c_i:
                inc = pt.elements.symbol(key).neutron.b_c_i
            else:
                inc = 0
                scattering_length += inc * 1j * component[key]
        scattering_length *= 1e-5
    return scattering_length
