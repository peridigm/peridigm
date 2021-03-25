from math import exp

class standard_linear_solid(object):
    def __init__(self, alpha, tau_b, modulus):
        self._alpha=alpha
        self._modulus=modulus
        self._tau_b=tau_b
        self._alpha_infinity=alpha * (1 - modulus)

    def get_alpha(self):
        return self._alpha

    def get_modulus(self):
        return self._modulus

    def get_tau_b(self):
        return self._tau_b

    def get_alpha_infinity(self):
        return self._alpha_infinity

    alpha=property(get_alpha)
    modulus=property(get_modulus)
    tau_b=property(get_tau_b)
    alpha_infinity=property(get_alpha_infinity)

    def relax(self, ed0, time_steps):
        td_long= ed0 * self.alpha_infinity
        analytical=lambda t: td_long + ed0 * self.modulus * self.alpha * exp(-t/self.tau_b)
        td=[analytical(t) for t in time_steps]
        return td
