from math import exp

class standard_linear_solid(object):
    def __init__(self, alpha, tau_b, tau):
        self._alpha=alpha
        self._tau=tau
        self._tau_b=tau_b
        self._e_infinity=alpha * (1 - tau_b/tau)
    
    def get_tau(self):
        return self._tau
    
    def get_tau_b(self):
        return self._tau_b
    
    def get_alpha(self):
        return self._alpha
    
    def get_e_infinity(self):
        return self._e_infinity
    
    alpha=property(get_alpha)
    tau=property(get_tau)
    tau_b=property(get_tau_b)
    e_infinity=property(get_e_infinity)
    
    def relax(self, ed0, time_steps):
        td_long= ed0 * self.e_infinity
        analytical=lambda t: td_long + self.alpha * self.tau_b * ed0 * exp(-t/self.tau_b) / self.tau 
        td=[analytical(t) for t in time_steps]
        return td
