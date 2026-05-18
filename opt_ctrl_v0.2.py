import numpy as np
import matplotlib.pyplot as plt
from track import track_segments, Track, Curve_BSpline
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, brute, Bounds
from models import m25_Final, m24
from plots import make_plots

model = m25_Final
data_file = './phase2/Gear4.csv'

L     = 0.15
m     = 0.18
gamma = 0.616 # 1.05
Cb    = 0.954 # 5.143
mu    = 1.75 # 0.355
g = 9.81
LAPS = 3

segments = track_segments['real']
trk = Track(segments)

def event_safety_violation(t, state, *args):
    q, q_dot, psi, psi_dot = state
   
    # Check lateral acceleration (measurable)
    kappa = trk.kappa(q) # self.curvatures[idx]
    lateral_accel = q_dot**2 * np.abs(kappa)
    muF = 1.75
    # Safety margin on physics limit
    return muF * g * 0.81 - lateral_accel

event_safety_violation.terminal = True
event_safety_violation.direction = -1

def find_best_gear_for_timeframe(i, last_gear, initial_state,
                                 event_horizont, t_frame):
    # Define the time span for the simulation (from t=0 to t=2)
    t_start = t_frame * i
    t_end   = t_start + (t_frame * event_horizont)


    t_span = (t_start, t_end)
    Nf = 50
    N = Nf * event_horizont + 1
    t_eval = np.linspace(t_start, t_end, int(N))

    # return 4, solve_ivp(
    #                 model, t_span, initial_state,
    #                 args = (trk, L, m, Fm_at(4),gamma, Cb, mu),
    #                 method='RK45',
    #                 events=[event_safety_violation],
    #                 t_eval=t_eval,
    #                 dense_output=True
    #             ), Nf
    best_gear = last_gear

    while True:
        sol = solve_ivp(
            model, t_span, initial_state,
            args = (trk, L, m, Fm_at(best_gear),gamma, Cb, mu),
            method='RK45',
            events=[event_safety_violation],
            t_eval=t_eval
        )

        if not sol.success or sol.t_events[0].size > 0:
            
            best_gear = best_gear - 1
            
            if best_gear == 0:
                print(f'Not solution found at segment from {t_start} to {t_end} seconds, for gear #{best_gear} using a timeframe of {t_frame} s with an horizon of {event_horizont}.')
                best_gear = last_gear
            
            if best_gear >= last_gear:
                return best_gear, solve_ivp(
                    model, t_span, initial_state,
                    args = (trk, L, m, Fm_at(best_gear),gamma, Cb, mu),
                    method='RK45',
                    events=[event_safety_violation],
                    t_eval=t_eval
                ), Nf
            del sol
            continue
            
        # if sol.success:
        if best_gear < last_gear or best_gear >= 8:
            return best_gear, sol, Nf
        if best_gear >= last_gear + 2:
            return best_gear, sol, Nf
        best_gear = best_gear + 1


def Fm_at(gear):
    a = 0.59142912
    b = 0.3773076
    calc_Fm = np.log(gear) * a + b
    # print(f'{calc_Fm}')
    return calc_Fm
         
def find_optimal(gear_init = 4, event_horizont = 1.5, t_frame = 0.25):
    tot_q  = np.array([])
    tot_qd = np.array([])
    tot_p  = np.array([])
    tot_pd = np.array([])
    tot_t  = np.array([])
    gear_timeframes = np.array([], dtype='uint8')
    gear = gear_init

    initial_q       = 0.0 # Initial position
    initial_q_dot   = 0.0 # Initial velocity  
    initial_psi     = 0.0 # Initial orientation angle
    initial_psi_dot = 0.0 # Initial angular velocity

    i = 0
    while True:
        initial_state = [initial_q, initial_q_dot, initial_psi, initial_psi_dot]

        # simulation = {
        #     't_start': t_start,
        #     't_end': t_end,
        #     'run': lambda gear: solve_ivp(
        #         model, t_span, initial_state,
        #         args = (
        #             trk, L, m, Fm_at(gear),
        #             gamma, Cb, mu
        #         ),
        #         method='RK45',
        #         events=[event_safety_violation],
        #         t_eval=t_eval
        #     ),
        # }

        gear, fit_sol, Nf = find_best_gear_for_timeframe(i, gear, initial_state, event_horizont, t_frame)

        if Nf >= len(fit_sol.y[0]):
            Nf = -1

        if fit_sol.t_events[0].size > 0: 
            return None
        tot_q  = np.append(tot_q,  fit_sol.y[0][:Nf])
        tot_qd = np.append(tot_qd, fit_sol.y[1][:Nf])
        tot_p  = np.append(tot_p,  fit_sol.y[2][:Nf])
        tot_pd = np.append(tot_pd, fit_sol.y[3][:Nf])
        tot_t  = np.append(tot_t,  fit_sol.t[:Nf])

        gear_timeframes = np.append(gear_timeframes, gear)

        # update for next iter
        initial_q       = fit_sol.y[0][Nf] # Initial position
        initial_q_dot   = fit_sol.y[1][Nf] # Initial velocity  
        initial_psi     = fit_sol.y[2][Nf] # Initial orientation angle
        initial_psi_dot = fit_sol.y[3][Nf] # Initial angular velocity
        i += 1

        t_start = t_frame * i
        t_end   = t_start + (t_frame)
        # print(f'from {t_start:2.3f} to {t_end:2.3f} - gear {gear}')

        if tot_q[-1] > LAPS * trk.length:
            break

    tot_q  = np.append(tot_q,  fit_sol.y[0][Nf+1])
    tot_qd = np.append(tot_qd, fit_sol.y[1][Nf+1])
    tot_p  = np.append(tot_p,  fit_sol.y[2][Nf+1])
    tot_pd = np.append(tot_pd, fit_sol.y[3][Nf+1])
    tot_t  = np.append(tot_t,  fit_sol.t[Nf+1])

    end_idx = np.searchsorted(tot_q, LAPS * trk.length)
    tot_q  = tot_q[:end_idx]
    tot_qd = tot_qd[:end_idx]
    tot_p  = tot_p[:end_idx]
    tot_pd = tot_pd[:end_idx]
    tot_t  = tot_t[:end_idx]


    phis = trk.phi(tot_q)
    cont_q = np.linspace(0,tot_q[-1],5000)
    kappas = trk.kappa(cont_q)

    muF = 1.75
    
    # vlim1 = np.sqrt(1 * 9.81 * muF * 0.26)
    # vlim18 = np.sqrt(0.8 * 9.81 * muF * 0.26)
    # vlim2 = np.sqrt(1 * 9.81 * muF * 0.34)
    # vlim28 = np.sqrt(0.8 * 9.81 * muF * 0.34)

    print(f'{gear_init}, {event_horizont}, {t_frame:0.2f}, {tot_t[-1]:0.3f}, {tot_q[-1]:0.3f}')
    print((gear, event_horizont, t_frame))
    print(f'sim_t={tot_t[-1]}, sim_q={tot_q[-1]}')
    vlim_dyn = np.sqrt(9.81 * muF / kappas)

    import matplotlib.pyplot as plt
    fig3 = plt.figure(figsize=(18, 8))
    plt.plot(cont_q, vlim_dyn, 
                color='red', linestyle='--', linewidth=2.5, label=f'Physics limit')
    plt.plot(cont_q, vlim_dyn * np.sqrt(0.81), 
                color='gray', linestyle='--', linewidth=2, label=f'Safety margin')

    # plt.plot(tot_q, tot_p, label='sim psi', color='blue')
    # plt.plot(tot_q, phis, label='trk phi', color='green')
    # plt.plot(tot_q, phis - tot_p, label='beta', color='darkviolet')
    # plt.plot(tot_q, tot_qd, label='sim_qd', color='orange')
    plt.plot(tot_q, tot_qd, linewidth=2.5, label=f'{t_frame} seconds', alpha=0.8)
    

    # plt.plot([tot_q[0],tot_q[-1]], [vlim1,vlim1], label='min_lim', color='darkviolet', ls='dashed')
    # plt.plot([tot_q[0],tot_q[-1]], [vlim18,vlim18], label='min_lim', color='darkviolet', ls='dashed')
    # plt.plot([tot_q[0],tot_q[-1]], [vlim2,vlim2], label='max_lim', color='blue', ls='dashed')
    # plt.plot([tot_q[0],tot_q[-1]], [vlim28,vlim28], label='max_lim', color='blue', ls='dashed')
    plt.legend()
    plt.ylim(0,3)
    plt.show()

    sol_m25 = own_sol(y=(tot_q, tot_qd, tot_p, tot_pd), t=tot_t, gears=gear_timeframes)
    sol_m24 = sol_m25
    config = {
        'save': False,
        'show': True,
        'zoom': False,
        'show data': True,
        'show m24': False,
        # 'comment': '',
        'print params': True,
    }
    
    # plt, ani = make_plots(segments, sol_m24, sol_m25, sol_m25.data_file, sol_m25.params, config)
    
    # if config['show']:
        # plt.show()
    
    return sol_m25


class own_sol():
    def __init__(self, y, t, gears):
        self.y = y
        self.t = t
        self.data_file = './phase2/Gear4.csv'
        self.params = (L, m, Fm_at(4),gamma, Cb, mu), (L, m, Fm_at(4),gamma, Cb, mu)
        self.gears = gears
        

# for c_gear in np.linspace(1, 5, 5):
#     for c_event_horizont in (3, 3.5, 4):#np.linspace(1.5, 2.5, 5):
#         for c_t_frame in (0.1, 0.25, 0.5, 1.0): #np.linspace(0.1, 1, 6):
#             try:        
#                sol = find_optimal(c_gear, c_event_horizont, c_t_frame)
#             except:
#                 pass


# def ocp(c_event_horizont):#, *args):
#     if type(c_event_horizont) == np.ndarray:
#         c_event_horizont = c_event_horizont.item()
#     sol = find_optimal(3, c_event_horizont, 0.02)
#     if sol == None:
#         return 1000
#     return sol.t[-1]
       

# ranges = (slice(1.5, 3, 0.05),)
# res = brute(ocp, ranges)

param_bounds = Bounds(
    [5*1.5,], # Lower bounds: Fm, gamma, Cb, mu
    [5*3.0,]  # Upper bounds: Fm, gamma, Cb, mu
)

# Model 2024 parameter bounds??
# param_bounds = Bounds(
#     [1.59, 0.30, 1.22, 0.65, 0.610], # Lower bounds: Fm, gamma, k, mu, g2
#     [1.59, 0.35, 1.26, 0.665, 0.616]  # Upper bounds: Fm, gamma, k, mu, g2
# )
# p0 = Fm, gamma, k, mu, g2

# res = minimize(ocp, 5*2.7, method='Nelder-Mead', bounds=param_bounds)#, options={'maxiter':20})
# popt = res.x
# # popt = ocp(2.7)
# print(popt)



print('Done!')

def plotter():
    gear = 3
    event_horiz = 3.8375
    t_frame = (0.02,0.1)#, 0.25, 0.5, 1.0)

    solutions = []

    solutions.append(find_optimal(3, 13.8375, 0.02))
    solutions.append(find_optimal(3, 3.8375, 0.1))

    # for tf in t_frame:
        # solutions.append(find_optimal(gear, event_horiz, tf))

    plot_lap_time_vs_control_points(gear, event_horiz, t_frame, solutions)
    plot_discrete_control_levels(gear, event_horiz, t_frame, solutions)
    plot_control_force_mapping(gear, event_horiz, t_frame, solutions)
    plot_control_force_mapping2(gear, event_horiz, t_frame, solutions)
    plot_speed_profiles(gear, event_horiz, t_frame, solutions)



def plot_lap_time_vs_control_points(gear, event_horiz, t_frame, solutions):
    fig = plt.figure(figsize=(10, 8))
    colors = ['blue', 'peru', 'green', 'darkviolet']

    for sol, tf, color in zip(solutions, t_frame, colors):
        lap_times = sol.t[-1] / LAPS

        plt.plot(tf, lap_times, 'o-', c=color, linewidth=3, markersize=12)
        plt.xlabel('Lookup Timeframe', fontsize=14)
        plt.ylabel('Lap Time [s]', fontsize=14)
        txt = plt.text(0.12, 5.3, f'Gear: {gear}\nEvent Horizon: {event_horiz}', fontsize=16, fontweight='bold')
        txt.set_bbox(dict(facecolor='lightgray', alpha=0.5, edgecolor='lightgray'))
    plt.title('Lap Time vs Lookup Timeframe', fontsize=16, fontweight='bold', pad=20)
    plt.grid(True, alpha=0.3)
    plt.xticks(t_frame, fontsize=12)
    plt.yticks(fontsize=12)
    plt.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.10)
    plt.show()


def plot_discrete_control_levels(gear, event_horiz, t_frame, solutions):
    fig = plt.figure(figsize=(15, 9))
    colors = ['blue', 'peru', 'green', 'darkviolet']

    for sol, tf, color in zip(solutions, t_frame, colors):
        t_step = event_horiz * tf
        levels = np.append(sol.gears, sol.gears[-1])
        times = np.append(np.linspace(0, len(sol.gears) * tf, len(sol.gears), endpoint=False), np.float64(sol.t[-1]))

        times_ext = np.insert(times,[0,times.size],[0,times[-1]])
        levels_ext = np.insert(levels,[0,levels.size],[levels[-1],levels[0]])

    
        plt.step(times_ext, levels_ext, '-', c=color, linewidth=3, where='post', label=f'{tf} seconds')
        # plt.scatter(times, levels,  s=200, zorder=5,  linewidth=1) #c='darkviolet', edgecolors='darkred',

    plt.xlabel('Time [s]', fontsize=14, labelpad=10)
    plt.ylabel('Discrete Control Level', fontsize=14, labelpad=10)
    plt.title('Discrete Control Strategy', fontsize=16, fontweight='bold', pad=25)
    plt.legend(fontsize=12)
    plt.ylim([-0.5, 8.5])
    plt.yticks(range(0, 8), fontsize=10)
    plt.xticks(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.subplots_adjust(left=0.08, right=0.96, top=0.92, bottom=0.10)
    plt.show()


def plot_control_force_mapping(gear, event_horiz, t_frame, solutions):
    fig = plt.figure(figsize=(14, 8))
    
    colors = ['blue', 'peru', 'green', 'darkviolet']

    for sol, tf, color in zip(solutions, t_frame, colors):
        t_step = event_horiz * tf
        levels = np.append(sol.gears, sol.gears[-1])
        times = np.append(np.linspace(0, len(sol.gears) * tf, len(sol.gears), endpoint=False), np.float64(sol.t[-1]))

        plt.step(times, Fm_at(levels), '-', c=color, linewidth=3, where='post',  label=f'{tf} seconds')
        # plt.plot(times, Fm_at(levels), linewidth=2.5, label=f'{tf} seconds', alpha=0.8)
    
    plt.xlabel('Time [s]', fontsize=14, labelpad=10)
    plt.ylabel('Motor Force [N]', fontsize=14)
    plt.title('Discrete Control Force Profiles', fontsize=16, fontweight='bold', pad=20)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.subplots_adjust(left=0.10, right=0.95, top=0.92, bottom=0.10)
    plt.show()


def plot_control_force_mapping2(gear, event_horiz, t_frame, solutions):
    fig = plt.figure(figsize=(14, 8))
    
    colors = ['blue', 'peru', 'green', 'darkviolet']

    for sol, tf, color in zip(solutions, t_frame, colors):
        end_idx = np.searchsorted(sol.y[0], trk.length * 2)

        positions = sol.y[0][:end_idx+1]
        time_space = np.append(np.linspace(0, len(sol.gears) * tf, len(sol.gears), endpoint=False), np.float64(sol.t[-1]))
        times = np.searchsorted(time_space, sol.t[:end_idx+1])
        levels = np.append(sol.gears, sol.gears[-1])
        gears = levels[times]

        plt.step(positions, Fm_at(gears), '-', c=color, linewidth=3, where='post',  label=f'{tf} seconds')
        # plt.plot(positions, Fm_at(levels), linewidth=2.5, label=f'{tf} seconds', alpha=0.8)
    
    plt.xlabel('Track Position [m]', fontsize=14, labelpad=10)
    plt.ylabel('Motor Force [N]', fontsize=14)
    plt.title('Discrete Control Force Profiles', fontsize=16, fontweight='bold', pad=20)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.subplots_adjust(left=0.10, right=0.95, top=0.92, bottom=0.10)
    plt.show()


def plot_speed_profiles(gear, event_horiz, t_frame, solutions):
    fig = plt.figure(figsize=(14, 8))

    colors = ['blue', 'peru', 'green', 'darkviolet']

    for sol, tf, color in zip(solutions, t_frame, colors):
        t_step = event_horiz * tf
        levels = np.append(sol.gears, sol.gears[-1])
        times = np.append(np.linspace(0, len(sol.gears) * tf, len(sol.gears), endpoint=False), np.float64(sol.t[-1]))
        # end_idx = np.searchsorted(times, sol.t[-1]) + 1

        # plt.figure()

        plt.plot(sol.y[0], sol.y[1], c=color, linewidth=2.5, label=f'{tf} seconds', alpha=0.8)
    
        cont_q = np.linspace(0, trk.length * LAPS, 5000)
        kappas = trk.kappa(cont_q)
        muF = 1.75
        vlim_dyn = np.sqrt(9.81 * muF / kappas)
    plt.plot(cont_q, vlim_dyn, 
                color='red', linestyle='--', linewidth=2.5, label=f'Physics limit')
    plt.plot(cont_q, vlim_dyn * np.sqrt(0.81), 
                color='gray', linestyle='--', linewidth=2, label=f'Safety margin')
           
    plt.xlabel('Track Position [m]', fontsize=14)
    plt.ylabel('Speed [m/s] (MEASURABLE)', fontsize=14)
    plt.title('Speed Profiles with Safety Validation\nUsing Only Measurable Quantities', fontsize=16, fontweight='bold', pad=20)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.xlim([0, trk.length])
    plt.ylim([1, 2.5])
    plt.show()


plotter()
find_optimal()