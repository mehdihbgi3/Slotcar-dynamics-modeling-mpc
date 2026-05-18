import numpy as np

def m24(time, state, trk, L, m, Fm, gamma, k, mu, gamma2):
    q, q_dot, psi, psi_dot=state
    
    if time > 3.14:
        pass
    dx=trk.dx(q) 
    dy=trk.dy(q) 
    ddx=trk.ddx(q) 
    ddy=trk.ddy(q) 
    phi=trk.phi(q)
    
    A=gamma*np.abs(np.cos(psi-phi))+gamma2*np.abs(np.sin(psi-phi))

    B=2.0*(dx**2+dy**2)
    C=L*(dx*np.sin(psi)-dy*np.cos(psi))
    D=2/m*Fm*(dx*np.cos(psi)+dy*np.sin(psi))-A*q_dot*(dx**2+dy**2)*4/m -2*q_dot**2*(dx*ddx+dy*ddy)-L*psi_dot**2*(dx*np.cos(psi)+dy*np.sin(psi))
    F=2/m*k*L*(phi-psi) -mu*2/m*L**2*psi_dot  -q_dot**2*L*(ddx*np.sin(psi)-ddy*np.cos(psi))
    E=L*(dx*np.sin(psi)-dy*np.cos(psi))

    q_double_dot=(D*L**2-C*F)/(B*L**2-C*E)
    psi_double_dot=(B*F-E*D)/(B*L**2-C*E)
    
    # Return the rates of change
    return np.array([q_dot, q_double_dot, psi_dot, psi_double_dot])

def m24_bspl(time, state, trk, L, m, Fm, gamma, k, mu, gamma2):
    q, q_dot, psi, psi_dot=state
    
    dx=trk.dx(q) 
    dy=trk.dy(q) 
    ddx=trk.ddx(q) 
    ddy=trk.ddy(q) 
    phi=trk.phi(q)
    
    A=gamma*np.abs(np.cos(psi-phi))+gamma2*np.abs(np.sin(psi-phi))

    B=2.0*(dx**2+dy**2)
    C=L*(dx*np.sin(psi)-dy*np.cos(psi))
    D=2/m*Fm*(dx*np.cos(psi)+dy*np.sin(psi))-A*q_dot*(dx**2+dy**2)*4/m -2*q_dot**2*(dx*ddx+dy*ddy)-L*psi_dot**2*(dx*np.cos(psi)+dy*np.sin(psi))
    F=2/m*k*L*(phi-psi) -mu*2/m*L**2*psi_dot-q_dot**2*L*(ddx*np.sin(psi)-ddy*np.cos(psi))
    E=L*(dx*np.sin(psi)-dy*np.cos(psi))

    q_double_dot=(D*L**2-C*F)/(B*L**2-C*E)
    psi_double_dot=(B*F-E*D)/(B*L**2-C*E)
    
    # Return the rates of change
    return np.array([q_dot, q_double_dot, psi_dot, psi_double_dot])

def m24_v2(time, state, trk, L, m, Fm, gamma, k, mu, gamma2):
    q, q_dot, psi, psi_dot = state
    
    dx = trk.dx(q) 
    dy = trk.dy(q) 
    ddx = trk.ddx(q) 
    ddy = trk.ddy(q) 
    phi = trk.phi(q)
    
    A = gamma * np.abs(np.cos(psi - phi)) + mu * np.abs(np.sin(psi - phi))

    B = 2.0 * (dx**2 + dy**2)
    C = L * (dx * np.sin(psi) - dy * np.cos(psi))
    E = C

    thrust_q = (Fm / m) * (dx * np.cos(psi) + dy * np.sin(psi))
    friction_q = -(A / m) * q_dot * (dx**2 + dy**2)
    inertial_q = -(2.0 * q_dot**2) * (dx * ddx + dy * ddy) \
                    - L * psi_dot**2 * (dx * np.cos(psi) + dy * np.sin(psi))

    D = 2.0 * thrust_q + 4.0 * friction_q + inertial_q

    spring_psi = (k / m) * L * (phi - psi)
    damping_psi = -(2.0 * mu / m) * (L**2) * psi_dot
    centrifugal_psi = -q_dot**2 * L * (ddx * np.sin(psi) - ddy * np.cos(psi))
    
    F = 2.0 * spring_psi + damping_psi + centrifugal_psi
    
    delta = B * L**2 - C * E
    q_double_dot = (D * L**2 - C * F) / delta
    psi_double_dot = (B * F - E * D) / delta
    
    # Return the rates of change
    return np.array([q_dot, q_double_dot, psi_dot, psi_double_dot])

def m25_Final(time, state, trk, L, m, Fm, gamma, Cb, Cr_mu):
    q, q_dot, psi, psi_dot = state
    
    dx = trk.dx(q) 
    dy = trk.dy(q) 
    ddx = trk.ddx(q) 
    ddy = trk.ddy(q) 
    phi = trk.phi(q)
    phi = np.mod(phi + np.pi, 2 * np.pi) - np.pi
    psi = np.mod(psi + np.pi, 2 * np.pi) - np.pi

    beta = phi - psi  # Slip angle
    beta = np.mod(beta + np.pi, 2 * np.pi) - np.pi

    B = 2.0 * (dx**2 + dy**2)
    C = L * (dx * np.sin(psi) - dy * np.cos(psi))

    # Motor thrust
    thrust_q = (Fm / m) * (dx * np.cos(psi) + dy * np.sin(psi))
    
    # Longitudinal friction (always opposes motion, never negative)
    # A = gamma - Cb * beta**2

    # A = gamma * np.abs(np.cos(beta)) + Cb * beta * np.sin(beta)
    A = gamma - Cb * beta **2
    if A < 0.0:
        A = 0.0
    
    friction_q = -(A / m) * q_dot * (dx**2 + dy**2) # t^2
    # Inertial terms
    inertial_q = -(2.0 * q_dot**2) * (dx * ddx + dy * ddy)
    centrifugal_psi = -L * psi_dot**2 * (dx * np.cos(psi) + dy * np.sin(psi))

    # D = 2.0 * thrust_q + 4.0 * friction_q + inertial_q
    Q_q = inertial_q + centrifugal_psi + thrust_q + friction_q 

    # Lateral cornering force (separate from friction)
    lateral_cornering = (Cb / m) * L * beta
    
    # Yaw damping
    damping_psi = -(Cr_mu / m) * (L**2) * psi_dot
    
    # Geometric coupling
    centrifugal_q = -q_dot**2 * L * (ddx * np.sin(psi) - ddy * np.cos(psi))

    Q_psi = centrifugal_q + lateral_cornering + damping_psi
    
    # Solve for accelerations
    delta = B * L**2 - C**2
    if abs(delta) < 1e-10:
        delta = np.sign(delta) * 1e-10
    
    q_double_dot = (Q_q * L**2 - C * Q_psi) / delta
    psi_double_dot = (B * Q_psi - C * Q_q) / delta
    
    # Clip for numerical stability
    q_double_dot = np.clip(q_double_dot, -30.0, 30.0)
    psi_double_dot = np.clip(psi_double_dot, -50.0, 50.0)
    
    # Return the rates of change
    return np.array([q_dot, q_double_dot, psi_dot, psi_double_dot])
