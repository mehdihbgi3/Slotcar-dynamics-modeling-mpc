import numpy as np
import matplotlib.pyplot as plt
from track import track_segments, Track, Curve_BSpline
from scipy.integrate import solve_ivp
from scipy.optimize import minimize, Bounds
from models import m25_Final, m24
from matplotlib.font_manager import FontProperties

model = m25_Final
data_file = './Gear5.csv'

L = 0.15
m = 0.18
car_param = L, m

exp_data = np.genfromtxt(data_file, delimiter=',')

exp_x   = exp_data[:,0]
exp_y   = exp_data[:,1]
exp_t   = exp_data[:,2]
exp_q   = exp_data[:,3]
exp_psi = exp_data[:,4]
exp_dq  = exp_data[:,5]

# Define the time span for the simulation (from t=0 to t=2)
t_span = (0.0, exp_t[-1])
print(t_span)
N = 5000
t_eval = np.linspace(t_span[0], t_span[1], N)

initial_q = 0.0         # Initial position
initial_q_dot = 0.0     # Initial velocity  
initial_psi = 0.0       # Initial orientation angle
initial_psi_dot = 0.0   # Initial angular velocity

y0 = [initial_q, initial_q_dot, initial_psi, initial_psi_dot]


# idx = [
#     (   0, 89),
#     ( 144, 310),
#     ( 339, 464),
#     ( 469, 490),
#     ( 501, 615),

# ]

# idx_start, idx_end = idx[1]
# # Define the time span for the simulation (from t=0 to t=2)
# t_start = exp_t[idx_start] # - t_offset
# t_end   = exp_t[idx_end] - t_start# - t_offset

# t_span = (0, t_end)
# N = int(100 * t_end)
# t_eval = np.linspace(t_span[0], t_span[1], N)

# initial_q = float(exp_q[idx_start])      # Initial position
# initial_q_dot = 0.0                      # Initial velocity  
# initial_psi = float(exp_psi[idx_start])  # Initial orientation angle
# initial_psi_dot = 0.0                    # Initial angular velocity

# y0 = [initial_q, initial_q_dot, initial_psi, initial_psi_dot]

segments = track_segments['real']
coupled_system = Track(segments)

def get_data(solution, trk):

    sim_q = np.interp(exp_t, solution.t, solution.y[0])
    sim_x = np.interp(sim_q % trk.length, trk.pt_q, trk.pt_x)
    sim_y = np.interp(sim_q % trk.length, trk.pt_q, trk.pt_y)
    sim_psi = np.interp(exp_t, solution.t, solution.y[2])
    sim_qd = np.interp(exp_t, solution.t, solution.y[1])

    return sim_x, sim_y, sim_q, sim_psi, sim_qd

def make_plot(params2):
    fit_sol = solve_ivp(model, t_span, y0, args=(coupled_system,*car_param,*params2), method='RK45', t_eval=t_eval)
    fm,gamma,Cb,mu = params2
    # fm,gamma,Cb,mu, gamma2 = params2
    

    sim_x, sim_y, sim_q, sim_psi, sim_qd = get_data(fit_sol, coupled_system)

    print(f'{sim_q[-1]=}, {exp_q[-1]=}, {exp_q[-1]/sim_q[-1]}')

    resid = np.sum((sim_x - exp_x)**2) + np.sum((sim_y - exp_y)**2)
    idxs = np.int64(np.linspace(0,len(sim_x)-1,9))
    dis = np.sum(np.hypot(sim_x[idxs]-exp_x[idxs], sim_y[idxs]-exp_y[idxs])).item()
    name = f'_{resid:.3f}_{dis=:.4f}_{mu=:.4f}_{gamma=:.4f}_{fm=:.4f}_{Cb=:.4f}.png'.replace('=','-')

    # if resid < 200:
    if True:
        dirpath = './plts/g1/'
    # elif resid < 300:
    #     dirpath = './plts/300/'
    # elif resid < 350:
    #     dirpath = './plts/350/'
    # elif resid < 400:
    #     dirpath = './plts/400/'
    else:
        return

    fig2 = plt.figure(figsize=(18, 8))
    plt.plot(exp_t, sim_x, label='sim X', color='blue')
    plt.plot(exp_t, exp_x, label='exp X', color='green')
    plt.plot(exp_t, sim_y, label='sim Y', color='red')
    plt.plot(exp_t, exp_y, label='exp Y', color='orange')

    plt.xlabel(r'$t$')
    plt.ylabel(r'$Z$')
    plt.grid()
    plt.legend(loc='upper left')
    plt.tight_layout()
    # plt.show()
    # plt.savefig(dirpath + 'plt' + name)
    plt.close(fig2)


    newdd = exp_q - sim_q % coupled_system.length
    newdd[newdd > 3] = newdd[newdd > 3]-coupled_system.length
    newdd[newdd < -3] = newdd[newdd < -3] + coupled_system.length

    fig1 = plt.figure(figsize=(18, 8))
    # plt.plot(exp_t,newdd)
    plt.plot(exp_t, sim_q % coupled_system.length, label='sim q', color='blue')
    plt.plot(exp_t, exp_q, label='exp q', color='green')

    plt.xlabel(r'$t$')
    plt.ylabel(r'$Z$')
    plt.grid()
    plt.legend(loc='upper left')
    plt.tight_layout()
    plt.show()
    # plt.savefig(dirpath + 'pltQ' + name)
    plt.close(fig1)



    # fig4 = plt.figure(figsize=(11.3, 3.2))
    # # fig3 = plt.figure(figsize=(7, 2.5))
    # plt.title(r"Heading angle $\psi$ vs time")
    # plt.xlabel("Time $t [s]$")
    # plt.ylabel(r"Heading angle $\psi [rad]$")

    # plt.plot(exp_t, exp_psi, color='red', label=r'Reference track angle $\phi$')
    # plt.plot(fit_sol.t, fit_sol.y[2], color='limegreen', label=r'Model 25 $\psi$')
    # plt.tight_layout()
    # plt.legend()
    # plt.savefig(dirpath + 'beta/' + name + '.png')
    # plt.close(fig4)




    fig3 = plt.figure(figsize=(18, 8))

    phis = coupled_system.phi(sim_q)

    muF = 1.75
    vlim1 = np.sqrt(1 * 9.81 * muF * 0.26)
    vlim18 = np.sqrt(0.8 * 9.81 * muF * 0.26)
    vlim2 = np.sqrt(1 * 9.81 * muF * 0.34)
    vlim28 = np.sqrt(0.8 * 9.81 * muF * 0.34)
    
    print(f'{fit_sol.t[-1]=}')

    # plt.plot(sim_q, sim_psi, label='sim psi', color='blue')
    # plt.plot(sim_q, phis, label='trk phi', color='green')
    # plt.plot(sim_q, phis - sim_psi, label='beta', color='red')
    # plt.plot(sim_q, sim_qd, label='sim_qd', color='orange')
    # plt.plot([sim_q[0],sim_q[-1]], [vlim1,vlim1], label='min_lim', color='red')
    # plt.plot([sim_q[0],sim_q[-1]], [vlim18,vlim18], label='min_lim', color='red')
    # plt.plot([sim_q[0],sim_q[-1]], [vlim2,vlim2], label='max_lim', color='blue')
    # plt.plot([sim_q[0],sim_q[-1]], [vlim28,vlim28], label='max_lim', color='blue')

    # plt.xlabel(r'$t$')
    # # plt.ylabel(r'$Z$')
    # plt.grid()
    # plt.legend(loc='upper left')
    # plt.tight_layout()
    # plt.xlim(1.5, 8)
    # # plt.ylim(1.5, 2.0)
    # plt.show()
    # # plt.savefig(dirpath + 'plt' + name)
    # plt.close(fig3)

gear = 5

def Fm_at(gear):
    a = 0.59142912
    b = 0.3773076
    calc_Fm = np.log(gear) * a + b
    print(f'{calc_Fm}')
    return calc_Fm

L     = 0.15
m     = 0.18
gamma = 0.616 # 1.05
Cb    = 0.954 # 5.143
mu    = 1.75#0.245 # 0.355
g = 9.81
LAPS = 1
# make_plot((0.388*1.21, gamma*1, Cb*0.25, mu))
# make_plot((0.65 * Fm_at(gear), 0.25*gamma, 0.5*Cb, 0.5*mu, 
#            0.25 * mu))
# make_plot((0.65 * Fm_at(gear), 0.154, 0.477, 0.875, 
#            0.25 * mu))
# for fm2 in np.arange(3*0.38, 3*0.38, 0.1):
    # make_plot((fm2, gamma, Cb, mu))
fm24 = 0.388
# for fm2 in np.linspace(0.6 * fm24, 0.65 * fm24, 11):
# for fm2 in np.arange(0.843, 0.845, 0.001):
#     for gamma in np.linspace(0.154,0.653,1):
#         for cb2 in np.linspace(0.477, 0.8,1):
#             for mu2 in np.linspace(0.375, 0.2,1):#np.linspace(0.426,0.426,1):#1.2):
#             # make_plot((float(gamma),float(mu2),float(fm2)))
#                 make_plot((
#                     float(fm2),
#                     float(gamma),
#                     float(cb2),
#                     0.875,
#                     float(mu2),
#                     ))

print('Done!')


def plot_speed_profiles(gears, colors):

    solutions = []

    # solutions.append(find_optimal(3, 13.8375, 0.02))
    # solutions.append(find_optimal(3, 3.8375, 0.1))

    for gear in gears:
        solutions.append(solve_ivp(
            model, t_span, y0,
            args = (coupled_system, L, m, Fm_at(gear),gamma, Cb, mu),
            method='RK45',
            t_eval=t_eval
        ))

    

    for sol, gear, color in zip(solutions, gears, colors):
        plt.plot(sol.y[0], sol.y[1], c=color, linewidth=2.5, label=f'Gear #{gear}', alpha=0.8)
    
    
        cont_q = np.linspace(0, coupled_system.length * LAPS, 5000)
        kappas = coupled_system.kappa(cont_q)
        muF = 1.75
        vlim_dyn = np.sqrt(9.81 * muF / kappas)

        x_dsc = np.searchsorted(sol.y[0], 2.01)
        x_sc = np.searchsorted(sol.y[0], 2.38)
        x_dhc = np.searchsorted(sol.y[0], 3.84)
        x_hc = np.searchsorted(sol.y[0], 4.28)

        plt.annotate(
            f't={sol.t[x_dsc]:.2f}', xy=(2.01, 1.92+0.1),
            fontproperties=FontProperties(size=16, weight='bold'),
            bbox=dict(facecolor='lightgray', alpha=0.5, edgecolor='lightgray')
        )
        plt.annotate(
            f't={sol.t[x_sc]:.2f}', xy=(2.38, 1.72+0.1), 
            fontproperties=FontProperties(size=16, weight='bold'),
            bbox=dict(facecolor='lightgray', alpha=0.5, edgecolor='lightgray')
        )
        plt.annotate(
            f't={sol.t[x_dhc]:.2f}', xy=(3.84, 1.97+0.1), 
            fontproperties=FontProperties(size=16, weight='bold'),
            bbox=dict(facecolor='lightgray', alpha=0.5, edgecolor='lightgray')
        )
        plt.scatter([2.01, 2.38, 3.84, 4.28], [1.9, 1.72, 1.97, 2.18], c='red', s=200, zorder=5, edgecolors='darkred', linewidth=1)
        plt.annotate(
            f't={sol.t[x_hc]:.2f}', xy=(4.28, 2.18+0.1), 
            fontproperties=FontProperties(size=16, weight='bold'),
            bbox=dict(facecolor='lightgray', alpha=0.5, edgecolor='lightgray')
        )



    plt.plot(cont_q, vlim_dyn, 
                color='red', linestyle='--', linewidth=2.5, label=f'Physics limit')
    plt.plot(cont_q, vlim_dyn * np.sqrt(0.81), 
                color='gray', linestyle='--', linewidth=2, label=f'Safety margin')
    
    vlim1 = np.sqrt(1 * 9.81 * muF * 0.26)
    vlim18 = np.sqrt(0.8 * 9.81 * muF * 0.26)

    plt.plot([2.38,2.38], [0,5], color='gray', linestyle='-', linewidth=2, label=f'Soft crash position')
    plt.plot([4.28,4.28], [0,5], color='red', linestyle='-', linewidth=2.5, label=f'Hard crash position')

    # plt.axvspan(2.38-0.25, 2.38, color='gray', alpha=0.5, lw=0)
    # plt.fill_between([2.38,2.38], 1, 5, facecolor='green', alpha=.5)
    # plt.fill_between(t, -1, where=s < 0, facecolor='red', alpha=.5)

           
    plt.xlabel('Track Position [m]', fontsize=14)
    plt.ylabel('Speed [m/s] (MEASURABLE)', fontsize=14)
    plt.title('Speed Profiles with Safety Validation\nUsing Only Measurable Quantities', fontsize=16, fontweight='bold', pad=20)
    plt.legend(fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.xlim([0, coupled_system.length])
    plt.ylim([1, 2.5])
    plt.show()


# plot_speed_profiles(
#     [1,2,3,4],
#     ['blue', 'peru', 'green', 'darkviolet', 'black']
# )
plot_speed_profiles(
    [5],
    ['darkviolet', 'black']
)
# plot_speed_profiles(
#     [4,5],
#     ['darkviolet', 'black']
# )