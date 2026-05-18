import numpy as np
import types
from scipy.interpolate import BSpline, splrep

# --- Track Definitions ---
class TrackSegment:
    _r_int = 0.26
    _r_ext = 0.34
    _straight_length = 0.3464
    
    def __init__(self, segment_type, qty=1.0, friction=1.0, slope=0.0, heigth = 0.0):
        if segment_type == 'S':  # Straight
            self.length = qty * self._straight_length
            self.radius = np.inf
            self.curvature = 1.5
        elif segment_type == 'R':  # Right turn
            self.length = (qty/3) * np.pi * self._r_ext
            self.radius = -self._r_ext
            self.curvature = abs(1 / self.radius)
        elif segment_type == 'L':  # Left turn
            self.length = (qty/3) * np.pi * self._r_int
            self.radius = self._r_int
            self.curvature = abs(1 / self.radius)
        else:
            raise ValueError("Unknown segment type. Use 'S', 'L', or 'R'.")
        self.slope = slope
        self.friction = friction
        self.type = segment_type
        self.heigth = heigth

    def set_heigth(self, segment):
        self.heigth = segment.length * segment.slope

def set_heigth(seg):
    for i in range(1,len(seg)):
        seg[i].set_heigth(seg[i-1])
    

car_params = {
    'length': 0.15,
    'mass': 0.18,
}

track_segments = {
    # 'straight': [TrackSegment('S', 6)],
    'straight': [TrackSegment('S', 10*5.813953488374/2)],
    'circle': [TrackSegment('R', 6)],
    'circle2': [TrackSegment('R', 1), TrackSegment('S', 4),],
    'real': [
        TrackSegment('S', 1.5),
        TrackSegment('R', 2),
        TrackSegment('S', 2),
        TrackSegment('L', 2),
        TrackSegment('S', 1.34),
        TrackSegment('S', 2.66),
        TrackSegment('L', 3),
        TrackSegment('S', 3.5, slope=0.1),
        TrackSegment('S', 2.84, slope=-0.1),
        TrackSegment('S', 0.66, slope=-0.1),
        TrackSegment('R', 1),
        TrackSegment('S', 1),
        TrackSegment('R', 2),
        TrackSegment('S', 1.5),
    ],
    'last_year': [
        TrackSegment('S', 4),
        TrackSegment('R', 3),
        TrackSegment('S', 4),
        TrackSegment('R', 3),
    ],
}

class Track:
    def __init__(self, track_segments: list[TrackSegment]):
        track_x, track_y, track_q, track_phi, track_idx, track_kappa = self.track_coordinates(track_segments)
        # track_y[-20:][track_y[-20:] > 0] = 0.0
        # track_x[-20:][track_x[-20:] > 0] = 0.0
        self.pt_x = track_x
        self.pt_y = track_y
        self.pt_q = track_q
        self.pt_phi = track_phi
        self.pt_kappa = track_kappa
        self.segments_index = track_idx
        
        self.length = self.pt_q[-1]

        self.pt_dx = np.gradient(self.pt_x, self.pt_q)
        self.pt_dy = np.gradient(self.pt_y, self.pt_q)
        self.pt_ddx = np.gradient(self.pt_dx, self.pt_q)
        self.pt_ddy = np.gradient(self.pt_dy, self.pt_q)

        # self.pt_dx  = self.make_d1(self.pt_x, self.pt_q)
        # self.pt_ddx = self.make_d2(self.pt_x, self.pt_q)
        # self.pt_dy  = self.make_d1(self.pt_y, self.pt_q)
        # self.pt_ddy = self.make_d2(self.pt_y, self.pt_q)

        self.x = lambda q: np.interp(q % self.length, self.pt_q, self.pt_x)
        self.y = lambda q: np.interp(q % self.length, self.pt_q, self.pt_y)

        self.phi   = lambda q: np.interp(q % self.length, self.pt_q, self.pt_phi)
        self.dx    = lambda q: np.interp(q % self.length, self.pt_q, self.pt_dx)
        self.ddx   = lambda q: np.interp(q % self.length, self.pt_q, self.pt_ddx)
        self.dy    = lambda q: np.interp(q % self.length, self.pt_q, self.pt_dy)
        self.ddy   = lambda q: np.interp(q % self.length, self.pt_q, self.pt_ddy)
        self.kappa = lambda q: self.pt_kappa[np.searchsorted(self.pt_q, q % self.length, side='right')]
        
    # def make_d1(self, z, q):
    #     dz = np.diff(np.hstack([z[-2] - z[-1], z]))
    #     dq = np.diff(np.hstack([q[-2] - q[-1], q]))
    #     return dz/dq

    # def make_d2(self, z, q):
    #     dzdq = self.make_d1(z, q)
    #     ddzdq = np.diff(np.hstack([dzdq, dzdq[0]]))
    #     dq = np.diff(np.hstack([q, q[1] + q[-1]]))
    #     return ddzdq/dq

    def track_coordinates(self, segments):
        x, y, s, phi = [0], [0], [0],  [0]
        idx = []
        kappa = []
        
        for seg in segments:
            idx.append(len(x)-1)
            resolution = int(seg.length/0.012)
            lds = np.linspace(0, seg.length, resolution)
            if seg.type == 'S':
                x.extend(x[-1] + lds[1:] * np.cos(phi[-1]))
                y.extend(y[-1] + lds[1:] * np.sin(phi[-1]))
                s.extend(s[-1] + lds[1:])
                phi.extend([phi[-1]]*(resolution-1))
                # kappa.extend([seg.curvature]*(resolution-1))
            else:
                angles = lds[1:] / seg.radius
                new_phi = phi[-1] + angles
                new_x = x[-1] - seg.radius * (np.sin(phi[-1]) - np.sin(new_phi))
                new_y = y[-1] + seg.radius * (np.cos(phi[-1]) - np.cos(new_phi))
                
                x.extend(new_x)
                y.extend(new_y)
                s.extend(s[-1] + lds[1:])
                phi.extend(new_phi)
            kappa.extend([seg.curvature]*(resolution-1))

        idx.append(len(x))
        kappa.append(seg.curvature)
        print(s[-1])
        return np.array(x), np.array(y), np.array(s), np.array(phi), np.array(idx), np.array(kappa)

def find_current_segment(s, segments):
    total_length = sum(seg.length for seg in segments)
    s_mod = s % total_length

    # Find current segment
    seg_idx = 0
    s_remain = s_mod
    for seg in segments:
        if s_remain < seg.length:
            break
        s_remain -= seg.length
        seg_idx += 1
    return segments[seg_idx % len(segments)]


# --- Track Coordinates ---
def track_coordinates(segments):
    x, y, s, phi = [0], [0], [0], [0]
    idx = []
    
    for seg in segments:
        idx.append(len(x)-1)
        resolution = int(seg.length/0.012)
        lds = np.linspace(0, seg.length, resolution)
        if seg.type == 'S':
            x.extend(x[-1] + lds[1:] * np.cos(phi[-1]))
            y.extend(y[-1] + lds[1:] * np.sin(phi[-1]))
            s.extend(s[-1] + lds[1:])
            phi.extend([phi[-1]]*(resolution-1))
        else:
            angles = lds[1:] / seg.radius
            new_phi = phi[-1] + angles
            new_x = x[-1] - seg.radius * (np.sin(phi[-1]) - np.sin(new_phi))
            new_y = y[-1] + seg.radius * (np.cos(phi[-1]) - np.cos(new_phi))
            
            x.extend(new_x)
            y.extend(new_y)
            s.extend(s[-1] + lds[1:])
            phi.extend(new_phi)

    idx.append(len(x))
    print(s[-1])
    return np.array(x), np.array(y), np.array(s), np.array(phi), idx

# --- Path to Global XY ---
def path_to_xy(s, segments):
    total_length = sum(seg.length for seg in segments)
    s = s % total_length
    x, y = [0], [0]
    phi = 0
    for seg in segments:
        if s < seg.length:
            if seg.type == 'S':
                dx = s * np.cos(phi)
                dy = s * np.sin(phi)
                return x[-1] + dx, y[-1] + dy
            else:
                r = seg.radius
                angle = s / abs(r)
                cx = x[-1] - r * np.sin(phi)
                cy = y[-1] + r * np.cos(phi)
                new_phi = phi + angle * np.sign(r)
                new_x = cx + r * np.sin(new_phi)
                new_y = cy - r * np.cos(new_phi)
                return new_x, new_y
        else:
            if seg.type == 'S':
                dx = seg.length * np.cos(phi)
                dy = seg.length * np.sin(phi)
                x.append(x[-1] + dx)
                y.append(y[-1] + dy)
            else:
                r = seg.radius
                angle = seg.length / abs(r)
                cx = x[-1] - r * np.sin(phi)
                cy = y[-1] + r * np.cos(phi)
                new_phi = phi + angle * np.sign(r)
                new_x = cx + r * np.sin(new_phi)
                new_y = cy - r * np.cos(new_phi)
                x.append(new_x)
                y.append(new_y)
                phi += angle * np.sign(r)
        s -= seg.length
    return x[-1], y[-1]


from scipy.interpolate import BSpline, interp1d
from scipy.integrate import cumulative_trapezoid

class Curve_BSpline:
    # def __init__(self, degree, Radius, a, points_on_semicircle, points_on_line):
    def __init__(self, track_segments: list[TrackSegment]):
        track_x, track_y, track_q, track_phi, track_idx = track_coordinates(track_segments)
        
        knot_vect = np.linspace(0,track_q[-1],len(track_x))
        self.x = BSpline(knot_vect, track_x, 4)
        self.y = BSpline(knot_vect, track_y, 4)
        
        t_values = np.linspace(0, track_q[-1], 500)
        pt_dx = self.x(t_values, 1)
        pt_dy = self.y(t_values, 1)
        self.pt_dq = np.sqrt(pt_dx**2 + pt_dy**2)

        self.pt_q = cumulative_trapezoid(self.pt_dq, t_values, initial=0)

        self.length = self.pt_q[-1]
        t_values = np.linspace(0, self.length, 500)

        self.dx  = lambda q: self.x(q % self.length, 1)
        self.ddx = lambda q: self.x(q % self.length, 2)
        self.dy  = lambda q: self.y(q % self.length, 1)
        self.ddy = lambda q: self.y(q % self.length, 2)

        dx = self.dx(t_values)
        dy = self.dy(t_values)

        pt_phi = np.arctan2(dy, dx)
        self.pt_phi = np.mod(pt_phi + np.pi, 2 * np.pi) - np.pi
        
        self.pt_phi2 = track_phi
        self.pp = np.interp(t_values, track_q, self.pt_phi2)

    def phi(self, q):
        q3 = np.array([
            (q - .15 + self.length) % self.length,
            q, 
            (q + .15 + self.length) % self.length,
        ])
        phi = np.interp(q3 % self.length, self.pt_q, self.pt_phi)
        phi2 = np.interp(q % self.length, self.pt_q, self.pp)
        test = phi.mean() - phi2
        if abs(test) > 0.001:
            pphi = phi.mean()
            pass
        return phi2
        

    def track_coordinates(self, segments):
        x, y, s, phi = [0], [0], [0],  [0]
        idx = []
        
        for seg in segments:
            idx.append(len(x)-1)
            resolution = int(seg.length/0.012)
            lds = np.linspace(0, seg.length, resolution)
            if seg.type == 'S':
                x.extend(x[-1] + lds[1:] * np.cos(phi[-1]))
                y.extend(y[-1] + lds[1:] * np.sin(phi[-1]))
                s.extend(s[-1] + lds[1:])
                phi.extend([phi[-1]]*(resolution-1))
            else:
                angles = lds[1:] / seg.radius
                new_phi = phi[-1] + angles
                new_x = x[-1] - seg.radius * (np.sin(phi[-1]) - np.sin(new_phi))
                new_y = y[-1] + seg.radius * (np.cos(phi[-1]) - np.cos(new_phi))
                
                x.extend(new_x)
                y.extend(new_y)
                s.extend(s[-1] + lds[1:])
                phi.extend(new_phi)

        idx.append(len(x))
        print(s[-1])
        return np.array(x), np.array(y), np.array(s), np.array(phi), idx

    
