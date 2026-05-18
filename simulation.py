import numpy as np
from track import car_params, track_segments, Track, Curve_BSpline
from scipy.integrate import solve_ivp

from models import m24, m25_Final
from plots import make_plots

initial_q       = 0.0 #.344 * 2  # Initial position
initial_q_dot   = 0.0 # 0.30   # Initial velocity  
initial_psi     = 0.0 # 0  # Initial orientation angle
initial_psi_dot = 0.0 # -3   # Initial angular velocity

y0 = [initial_q, initial_q_dot, initial_psi, initial_psi_dot]

segments = track_segments['real']
# segments = track_segments['last_year']
# segments = track_segments['circle2']
# segments = track_segments['straight']
coupled_system = Track(segments)
# coupled_system = Curve_BSpline(segments)

L = car_params['length'] # 0.15
m = car_params['mass'] # 0.18




class Race:
    def __init__(self, y0, models, params, data_file, track):

        exp_data = np.genfromtxt(data_file, delimiter=',')
        exp_t   = exp_data[:,2]

        t_start = 0.0
        t_end = exp_t[-1]

        self.t_span = (t_start, t_end)
        self.t_eval = np.append(np.linspace(t_start, t_end, int(t_end/0.0025), endpoint=False), np.float64(t_end))
        self.y0 = y0
        self.m24 = models[0]
        self.m25 = models[1]
        self.params_24 = params[0]
        self.params_25 = params[1]
        self.data_file = data_file
        self.trk = track
    
    def solve(self, config):
        sol_m24 = solve_ivp(self.m24, self.t_span, self.y0, args=(self.trk, *self.params_24), method='RK45', t_eval=self.t_eval)
        sol_m25 = solve_ivp(self.m25, self.t_span, self.y0, args=(self.trk, *self.params_25), method='RK45', t_eval=self.t_eval)

        print(sol_m24.message)
        print(sol_m25.message)

        print(f"M24 q_end = {sol_m24.y[0][-1]}")
        print(f"M25 q_end = {sol_m25.y[0][-1]}")

        plt, ani = make_plots(segments, sol_m24, sol_m25, self.data_file, (self.params_24, self.params_25), config)
        if config['show']:
            plt.show()

        return


def params_m25(gear):
    # Model 25 parameters
    Fm_25 = [0.0, 0.388, 0.765, 1.036, 1.191, 1.338]
    gamma = 0.616
    Cb =    0.954
    mu =    1.75#0.245
    return L, m, Fm_25[gear], gamma, Cb, mu

def params_m24(gear, prepro='OK'):
    # Model 24 parameters
    Fm_24 = [0.0, 0.249, 0.494, 0.653, 0.765, 1.52]
    gamma_24 = 0.616/4.0
    k_24 = 0.954/2.0
    mu_24 = 1.75/2.0
    gamma2 = 0.375

    if prepro == 'OK':
        return L, m, Fm_24[gear], gamma_24, k_24, mu_24, gamma2
    
    Fm_24_w_bspl = [0.0, 0.0, 0.0, 1.315, 1.646, 1.52]
    gamma_24_w_bspl = 0.326
    k_24_w_bspl = 1.256
    mu_24_w_bspl = 0.65
    gamma2_w_bspl = 0.612

    # if prepro == 'bspl':
    #     return L, m, Fm_24_w_bspl[gear], gamma_24_w_bspl, k_24_w_bspl, mu_24_w_bspl, gamma2_w_bspl
    
    # # Model 24 parameters w/o data prepros
    # Fm_24_wo_dSPL = [0.0, 0.0, 0.0, 1.05, 1.364, 1.52]

    # Fm3_24_wo_dp = 1.0
    # Fm4_24_wo_dp = 1.311
    # Fm5_24_wo_dp = 1.52

    # gamma_24_wo_dp = 0.265
    # k_24_wo_dp = 1.27
    # mu_24_wo_dp = 0.685
    # gamma2_wo_dp = 0.61

    # if prepro == 'wo':
    #     return L, m, Fm_24_wo_dSPL[gear], gamma_24_wo_dp, k_24_wo_dp, mu_24_wo_dp, gamma2_wo_dp
    



config = {
    'save': False,
    'show': False,
    'zoom': False,
    'show data': False,
    'show m24': False,
    # 'comment': '',
    'print params': True,
}

gear1_race = Race(
    y0 = y0,
    models = (m24, m25_Final),
    params = (
        params_m24(1,'OK'), 
        params_m25(1)
    ),
    data_file = './Gear1.csv',
    track = coupled_system
)

gear2_race = Race(
    y0 = y0,
    models = (m24, m25_Final),
    params = (
        params_m24(2,'OK'), 
        params_m25(2)
    ),
    data_file = './Gear2.csv',
    track = coupled_system
)

gear3_race = Race(
    y0 = y0,
    models = (m24, m25_Final),
    params = (
        params_m24(3,'OK'), 
        params_m25(3)
    ),
    data_file = './Gear3.csv',
    track = coupled_system
)

gear4_race = Race(
    y0 = y0,
    models = (m24, m25_Final),
    params = (
        params_m24(4,'OK'), 
        params_m25(4)
    ),
    data_file = './Gear4.csv',
    track = coupled_system
)

gear5_race = Race(
    y0 = y0,
    models = (m24, m25_Final),
    params = (
        params_m24(5,'OK'), 
        params_m25(5)
    ),
    data_file = './phase1/data_set_old/gear5_full.csv',
    track = coupled_system
)

gear1_race.solve(config)
gear2_race.solve(config)
gear3_race.solve(config)
gear4_race.solve(config)

# I have not fitted gear 5 yet
# gear5_race.solve(config)
