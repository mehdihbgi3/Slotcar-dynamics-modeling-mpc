import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.widgets import Slider
from track import car_params, find_current_segment, Track

def get_kappas(qs, segments):
    segs = [find_current_segment(q, segments) for q in qs]
    return -np.array([seg.curvature for seg in segs])

def make_plots(segments, sol_m24, sol_m25, data_file, params, config):

    # comment = config['comment'] or ''

    
    if config['print params']:
        Fm, gam, Cb, Cr = params[1][2:]
        Fm /= 1.65
        gam /= 1.032
        Cb /= 5.045
        Cr /= 0.337 

        name = f'_{Fm=:.1f}_{gam=:.1f}_{Cb=:.1f}_{Cr=:.1f}.png'.replace('=','-') 
        # name = f'_{Fm=:.3f}_{gam=:.3f}_{Cb=:.3f}_{Cr=:.3f}.png'.replace('=','-') 
    else:
        name = ''
    # name = ''
    
    trk = Track(segments)
    track_x, track_y, track_s, track_phi, _, _ = trk.track_coordinates(segments)
    m25_q, m25_q_dot, m25_psi, m25_psi_dot = sol_m25.y
    m24_q, m24_q_dot, m24_psi, m24_psi_dot = sol_m24.y

    exp_data = np.genfromtxt(data_file, delimiter=',')
    exp_x   = exp_data[:,0]
    exp_y   = exp_data[:,1]
    exp_t   = exp_data[:,2]
    exp_q   = exp_data[:,3]
    exp_psi = exp_data[:,4]
    exp_dq  = exp_data[:,5]

    xx_exp = np.interp(sol_m25.t, exp_t, exp_x)
    yy_exp = np.interp(sol_m25.t, exp_t, exp_y)
    exp_qi = np.interp(sol_m25.t, exp_t, exp_q)
    
    xx_m25 = np.interp(m25_q % trk.length, trk.pt_q, trk.pt_x)
    yy_m25 = np.interp(m25_q % trk.length, trk.pt_q, trk.pt_y)
    xx_m24 = np.interp(m24_q % trk.length, trk.pt_q, trk.pt_x)
    yy_m24 = np.interp(m24_q % trk.length, trk.pt_q, trk.pt_y)
    
    # fig0 = plt.figure(figsize=(11.3, 3.2))
    # fig0 = plt.figure(figsize=(7, 2.5))
    # plt.title("Velocity vs time")
    # plt.xlabel("Time $t [s]$")
    # plt.ylabel("Vel $x|y [m]$")

    # # plt.plot(sol_m25.t, m25_q_dot, color='limegreen', label=r'$\dot{q} \text{ with } T_{tr,q}$')
    # # plt.plot(sol_m25.t, m25_q_dot, color='limegreen', label=r'$\dot{q} \text{ with } T_{rot,q}$')
    # plt.plot(sol_m25.t, m25_q_dot, color='limegreen', label=r'$\dot{q} \text{ with } T_{rot,\psi}$')
    # if config['show m24']:
    #     plt.plot(sol_m25.t, m24_q_dot, color='blue', label='Constant Speed')
    
    # plt.tight_layout()
    # plt.legend()
    # plt.savefig('./phase2/plts/T_rot_q__q_dot.png')
    # plt.savefig('./phase2/plts/T_rot_psi__q_dot.png')
    # plt.savefig('./phase2/plts/T_tr_q__q_dot.png')
    
    fig1 = plt.figure(figsize=(11.3, 3.2))
    plt.title("Position vs time")
    plt.xlabel("Time $t [s]$")
    plt.ylabel("Position $x|y [m]$")

    if config['show data']:
        plt.plot(sol_m25.t, xx_exp, color='red', label='Exp. data $x$')
        plt.plot(sol_m25.t, yy_exp, color='brown', label='Exp. data $y$')
    plt.plot(sol_m25.t, xx_m25, color='limegreen', label='Model 25 $x$')
    plt.plot(sol_m25.t, yy_m25, color='green', label='Model 25 $y$')
    if config['show m24']:
        plt.plot(sol_m25.t, xx_m24, color='blue', label='OCP $x$')
        plt.plot(sol_m25.t, yy_m24, color='dodgerblue', label='OCP $y$')
    
    plt.tight_layout()
    plt.legend()
    # plt.ylim(0.0, 2.5)

    if config['save']:
        plt.savefig('./phase2/plts/res_xy_1F'+name+'.png')
        if config['zoom']:
            plt.xlim(0,8)
            plt.savefig('./phase2/plts/res_xy_2Z'+name+'.png')
        plt.close(fig1)


    # fig2 = plt.figure(figsize=(11.3, 3.2))
    # plt.title("Position $Q$ vs time")
    # plt.xlabel("Time $t [s]$")
    # plt.ylabel("Position $q [m]$")

    # if config['show data']:
    #     plt.plot(sol_m25.t, exp_qi, color='red', label='Exp. data $q$')
    # plt.plot(sol_m25.t, m25_q, color='limegreen', label='Model 25 $q$')
    # if config['show m24']:
    #     plt.plot(sol_m25.t, m24_q, color='blue', label='OCP $q$')
    
    # plt.tight_layout()
    # plt.legend()

    # if config['save']:
    #     plt.savefig('./phase2/plts/res_q_1F'+name+'.png')
    #     if config['zoom']:
    #         plt.xlim(0,8)
    #         plt.savefig('./phase2/plts/res_q_2Z'+name+'.png')
    #     plt.close(fig2)
    

    # # fig3 = plt.figure(figsize=(11.3, 3.2))
    # fig3a = plt.figure(figsize=(7, 2.5))
    # plt.title(r"Angular velocity $\dot{\psi}$ vs time")
    # plt.xlabel("Time $t [s]$")
    # plt.ylabel(r"Angular velocity $\dot{\psi} [rad/s]$")

    # if config['show data']:
    #     plt.plot(sol_m25.t, trk.phi(exp_qi), color='red', label=r'Reference track angle $\phi$')

    # # plt.plot(sol_m25.t, m25_psi_dot, color='limegreen', label=r'$\dot{\psi} \text{ with } T_{tr,q}$')
    # # plt.plot(sol_m25.t, m25_psi_dot, color='limegreen', label=r'$\dot{\psi} \text{ with } T_{rot,q}$')
    # plt.plot(sol_m25.t, m25_psi_dot, color='limegreen', label=r'$\dot{\psi} \text{ with } T_{rot,\psi}$')
    # if config['show m24']:
    #     plt.plot(sol_m25.t, m24_psi_dot, color='blue', label=r'Constant Speed')

    # plt.tight_layout()
    # plt.legend()
    # # plt.savefig('./phase2/plts/centrifugal_psi.png')
    # # plt.savefig('./phase2/plts/T_rot_q__psi_dot.png')
    # # plt.savefig('./phase2/plts/T_rot_psi__psi_dot.png')
    # # plt.savefig('./phase2/plts/T_tr_q__psi_dot.png')


    # fig3 = plt.figure(figsize=(11.3, 3.2))
    # # fig3 = plt.figure(figsize=(7, 2.5))
    # plt.title(r"Heading angle $\psi$ vs time")
    # plt.xlabel("Time $t [s]$")
    # plt.ylabel(r"Heading angle $\psi [rad]$")

    # if config['show data']:
    #     plt.plot(sol_m25.t, trk.phi(exp_qi), color='red', label=r'Reference track angle $\phi$')

    # # plt.plot(sol_m25.t, m25_psi, color='limegreen', label=r'$\psi \text{ with } T_{tr,q}$')
    # # plt.plot(sol_m25.t, m25_psi, color='limegreen', label=r'$\psi \text{ with } T_{rot,q}$')
    # # plt.plot(sol_m25.t, m25_psi, color='limegreen', label=r'$\psi \text{ with } T_{rot,\psi}$')
    # plt.plot(sol_m25.t, m25_psi, color='limegreen', label=r'Model 25 $\psi$')
    # if config['show m24']:
    #     # plt.plot(sol_m25.t, m24_psi, color='blue', label=r'Constant Speed')
    #     plt.plot(sol_m25.t, m24_psi, color='blue', label=r'OCP $\psi$')

    # plt.tight_layout()
    # plt.legend()

    # # plt.savefig('./phase2/plts/T_rot_q__psi.png')
    # # plt.savefig('./phase2/plts/T_rot_psi__psi.png')
    # # plt.savefig('./phase2/plts/T_tr_q__psi.png')

    # # plt.ylim(-0.5, 3.5)

    # if config['save']:
    #     plt.savefig('./phase2/plts/res_psi_1F'+name+'.png')
    #     if config['zoom']:
    #         plt.xlim(0,8)
    #         plt.savefig('./phase2/plts/res_psi_2Z'+name+'.png')
    #     plt.close(fig3)

    # fig4 = plt.figure(figsize=(11.3, 3.2))
    # plt.title(r"Heading angle $\psi$ vs time - Residuals" + "\n" + r"i.e. Sideslip angle $\beta [rad]$")
    # plt.xlabel("Time $t [s]$")
    # plt.ylabel(r"Sideslip angle $\beta [rad]$")

    # plt.plot(sol_m25.t, m25_psi - trk.phi(exp_qi), color='limegreen', label=r'M25 $\beta$ Residuals') 
    # if config['show m24']:
    #     plt.plot(sol_m25.t, m24_psi - trk.phi(exp_qi), color='blue', label=r'M24 $\beta$ Residuals')    
    
    # plt.tight_layout()
    # plt.legend()

    # if config['save']:
    #     plt.savefig('./phase2/plts/res_beta_1F'+name+'.png')
    #     if config['zoom']:
    #         plt.xlim(0,8)
    #         plt.savefig('./phase2/plts/res_beta_2Z'+name+'.png')
    #     plt.close(fig4)

    fig5 = plt.figure(figsize=(11.3, 3.2))
    plt.title("Position $Q$ vs time - Residuals")
    plt.xlabel("Time $t [s]$")
    plt.ylabel("Position $q [m]$")
    
    m25_q_res = (m25_q % trk.length - exp_qi)
    m25_q_res[m25_q_res > 2] = m25_q_res[m25_q_res > 2] - trk.length
    m25_q_res[m25_q_res < -2] = m25_q_res[m25_q_res < -2] + trk.length
    m25_q_res2 = np.diff(m25_q_res, prepend=0)
    m25_q_res2[np.abs(m25_q_res2) > 0.15] = 0
    m25_q_res2 = np.cumsum(m25_q_res2)
    # m25_q_res2[m25_q_res2 > 6] = m25_q_res2[m25_q_res2 > 6] - trk.length
    # 
    # m25_q_res = np.convolve(m25_q_res, np.ones(5)/5, mode='valid')

    
    plt.plot(sol_m25.t,  m25_q_res2, color='limegreen', label='Model 25 $q$ Residuals')
    # plt.axhline(np.mean(m25_q_res2), color='limegreen', label='Model 25 $q$ Residuals')
    
    if config['show m24']:
        m24_q_res = (m24_q % trk.length - exp_qi) % trk.length
        m24_q_res[m24_q_res > 2] = m24_q_res[m24_q_res > 2] - trk.length
        m24_q_res[m24_q_res < -2] = m24_q_res[m24_q_res < -2] + trk.length
        m24_q_res2 = np.diff(m24_q_res, prepend=0)
        m24_q_res2[np.abs(m24_q_res2) > 0.15] = 0
        m24_q_res2 = np.cumsum(m24_q_res2)
        m24_q_res2[m24_q_res2 > 6] = m24_q_res2[m24_q_res2 > 6] - trk.length
        plt.plot(sol_m25.t, m24_q_res2,  color='blue', label='OCP $q$ Residuals')
        # plt.axhline(np.mean(m24_q_res2), color='blue', label='OCP $q$ Residuals')
    
    plt.tight_layout()
    plt.legend()

    if config['save']:
        plt.savefig('./phase2/plts/res_q_Res_1F'+name+'.png')
        if config['zoom']:
            plt.xlim(0,8)
            plt.savefig('./phase2/plts/res_q_Res_2Z'+name+'.png')
        plt.close(fig5)

    # print(f"Exp q_end = {exp_qi[-1]}")
    print(f"M24 RMSE = {np.sqrt(np.mean((m24_q - exp_qi)**2))}")
    print(f"M25 RMSE = {np.sqrt(np.mean((m25_q - exp_qi)**2))}")

    rmse_x_m24 = np.mean((xx_m24 - xx_exp)**2)
    rmse_y_m24 = np.mean((yy_m24 - yy_exp)**2)
    rmse_x_m25 = np.mean((xx_m25 - xx_exp)**2)
    rmse_y_m25 = np.mean((yy_m25 - yy_exp)**2)

    print(f"M24 RMSE xy = {np.sqrt(np.mean([rmse_x_m24, rmse_y_m24]))}")
    print(f"M25 RMSE xy = {np.sqrt(np.mean([rmse_x_m25, rmse_y_m25]))}")


    # fig6 = plt.figure(figsize=(11.3, 3.2))
    # plt.title("Position vs time - Residuals")
    # plt.xlabel("Time $t [s]$")
    # plt.ylabel("Position $x|y [m]$")

    # plt.plot(sol_m25.t, xx_m25 - xx_exp, color='limegreen', label='Model 25 $x$ Residuals')
    # plt.plot(sol_m25.t, yy_m25 - yy_exp, color='green', label='Model 25 $y$ Residuals')
    # if config['show m24']:
    #     plt.plot(sol_m25.t, xx_m24 - xx_exp, color='blue', label='OCP $x$ Residuals')
    #     plt.plot(sol_m25.t, yy_m24 - yy_exp, color='dodgerblue', label='OCP $y$ Residuals')

    # plt.tight_layout()
    # plt.legend()

    # if config['save']:
    #     plt.savefig('./phase2/plts/res_xy_Res_1F'+name+'.png')
    #     if config['zoom']:
    #         plt.xlim(0,8)
    #         plt.savefig('./phase2/plts/res_xy_Res_2Z'+name+'.png')
    #     plt.close(fig6)
    
    # --- Updated Animation for Two-Point Car Representation ---
    fig, axo= plt.subplots(1, figsize=(6.5, 3.5))

    axo.plot(track_x, track_y, 'k--', alpha=0.4)
    
    exp = {}
    if config['show data']:
        exp = {
            'Wheel_F': axo.plot([], [], color='red', marker='o', mec='k', markersize=6, label='Data')[0],
            'Wheel_R': axo.plot([], [], color='red', marker='o', mec='k', markersize=6)[0],
            'Car_Body': axo.plot([], [], color='red', ls='-', lw=3)[0]
        }

    m25 = {
        # 'Wheel_F': axo.plot([], [], color='limegreen', marker='o', mec='k', markersize=6, label=r'$\text{Effect of } T_{rot,q}$')[0],
        # 'Wheel_F': axo.plot([], [], color='limegreen', marker='o', mec='k', markersize=6, label=r'$\text{Effect of } T_{tr,q}$')[0],
        # 'Wheel_F': axo.plot([], [], color='limegreen', marker='o', mec='k', markersize=6, label=r'$\text{Effect of } T_{rot,\psi}$')[0],
        'Wheel_F': axo.plot([], [], color='limegreen', marker='o', mec='k', markersize=6, label='M25')[0],
        'Wheel_R': axo.plot([], [], color='limegreen', marker='o', mec='k', markersize=6)[0],
        'Car_Body': axo.plot([], [], color='limegreen', ls='-', lw=3)[0]
    }

    m24 = {
        # 'Wheel_F': axo.plot([], [], color='blue', marker='o', mec='k', markersize=6, label='Constant Velocity')[0],
        'Wheel_F': axo.plot([], [], color='blue', marker='o', mec='k', markersize=6, label='M24')[0],
        'Wheel_R': axo.plot([], [], color='blue', marker='o', mec='k', markersize=6)[0],
        'Car_Body': axo.plot([], [], color='blue', ls='-', lw=3)[0]
    }

    time_text  = axo.text(0.02, 0.7, '', transform=axo.transAxes)

    axo.set_aspect('equal')
    # axo.set_xlim(min(track_x) - 0.75, max(track_x) - 0.25)
    # axo.set_ylim(min(track_y) - 0.5, max(track_y) + 0.5)
    axo.set_xlim(min(track_x) - 0.5, max(track_x) + 0.5)
    axo.set_ylim(min(track_y) - 0.5, max(track_y) + 0.5)
    # axo.set_xlim(-1.2, 1.6)
    # axo.set_ylim(-1.0, 0.3)
    axo.set_title("Slot Car as Two Points: Front Fixed to Track")
    axo.legend()

    qs = [m25_q]
    psis = [m25_psi]
    models = [m25]

    if config['show data']:
        qs.insert(0, exp_qi)
        psis.insert(0, trk.phi(exp_qi))
        models.insert(0, exp)
    if config['show m24']:
        qs.append(m24_q)
        psis.append(m24_psi)
        models.append(m24)


    # qs = (exp_qi, m25_q, m24_q)
    # psis = (trk.phi(exp_qi), m25_psi, m24_psi)
    # models = (exp, m25, m24)
    L = car_params['length']

    def animate(i):
        for q, psi, model in zip(qs, psis, models):
            x_f = trk.x(q[i])# - 0.5
            y_f = trk.y(q[i])

            x_r = x_f - L * np.cos(psi[i])
            y_r = y_f - L * np.sin(psi[i])

            model['Wheel_F'].set_data([x_f], [y_f])
            model['Wheel_R'].set_data([x_r], [y_r])
            model['Car_Body'].set_data([x_r, x_f], [y_r, y_f])
        
        time_text.set_text(f't = {sol_m25.t[i]:.6f}')

        ret_list = []
        if config['show data']:
            ret_list.extend([exp['Wheel_F'], exp['Wheel_R'], exp['Car_Body']])
        
        ret_list.extend([m25['Wheel_F'], m25['Wheel_R'], m25['Car_Body']])

        if config['show m24']:
            ret_list.extend([m24['Wheel_F'], m24['Wheel_R'], m24['Car_Body']])

        ret_list.append(time_text)

        return ret_list

    ani = animation.FuncAnimation(fig, animate, frames=range(0,len(sol_m25.t),16), interval=40, blit=True, cache_frame_data=True)

    # ani.save('./phase2/plts/g4_ok.mp4', fps=25)
    # ani.save('./phase2/plts/centrifugal_v.mp4', fps=25)

    # ani.save('./phase2/plts/T_tr_q.mp4', fps=25)
    # ani.save('./phase2/plts/T_rot_q.mp4', fps=25)
    # ani.save('./phase2/plts/T_rot_psi.mp4', fps=25)

    return plt, ani

