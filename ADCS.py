from math import pi, sqrt, ceil
from main import *

class Satellite:
    def __init__(self, inertia, max_ampl_disturbance_torque, thrust_moment_arm, manoeuvre_degrees, manoeuvre_secs, pointing_knowledge, orbit_altitude, planet: str):
        self.inertia = inertia # kg m²
        self.max_ampl_disturbance_torque = max_ampl_disturbance_torque # N m
        self.thrust_moment_arm = thrust_moment_arm  # m
        self.manoeuvre_degrees = manoeuvre_degrees # deg
        self.manoeuvre_secs = manoeuvre_secs # s
        self.pointing_knowledge = pointing_knowledge # arcsec
        self.pointing_knowledge_degrees = pointing_knowledge/3600 # deg
        self.orbit_dist = orbit_altitude + planet_data[planet]['mean_radius'] # m
        self.orbital_period = 2*pi*sqrt(self.orbit_dist**3/planet_data[planet]['gravitational_parameter'])

    def get_min_h_for_gyro_stiffness(self):
        '''
        Returns the minimum momentum a wheel needs to have to neutralise periodic disturbance torques to within the allowed range
        '''
        self.min_h = self.max_ampl_disturbance_torque * self.orbital_period / (self.pointing_knowledge_degrees*pi/180 * 4) # kg m² / s
        self.min_omega = self.min_h / self.wheel_inertia
        return self.min_h
    
    def get_max_h_in_periodic_disturbance(self):
        '''
        Returns the maximum momentum a reaction wheel will have in counteracting the periodic disturbance torque
        '''
        self.h_disturb = sqrt(2)/2 * self.max_ampl_disturbance_torque * self.orbital_period/4 # kg m² / s
        self.omega_disturb = self.h_disturb / self.wheel_inertia
        return self.h_disturb
    
    def define_wheel(self, inertia, density, thickness_ratio, inner_outer_ratio):
        '''
        Calculate the properties of a wheel from its inertia, density, thickness (height/diameter) and ratio of inner radius to outer radius (if not solid)
        '''
        self.wheel_inertia = inertia # kg m²
        self.wheel_density = density # kg/m³
        self.wheel_r_outer = (inertia / ((1-inner_outer_ratio**4)*pi*thickness_ratio*density))**(1/5) # m
        self.wheel_r_inner = inner_outer_ratio * self.wheel_r_outer # m
        self.wheel_mass = 2 * pi * density * self.wheel_r_outer**3 * thickness_ratio * (1-inner_outer_ratio**2) # kg
        
    def get_slew_torque(self):
        '''
        Returns the torque needed to satisfy slew requirements
        '''
        self.slew_torque = 4 * self.inertia * self.manoeuvre_degrees / self.manoeuvre_secs**2 # N m
        return self.slew_torque

    def thrust_for_wheel_dump(self, h, time):
        '''
        Returns thrust to dump momentum h in 'time' seconds
        '''
        F = h/(time*self.thrust_moment_arm) # N
        return F
    
steel_density = 8050 # kg/m³

leo = Satellite(inertia=200, 
                max_ampl_disturbance_torque=1e-2,
                thrust_moment_arm=1,
                manoeuvre_degrees=25,
                manoeuvre_secs=60,
                pointing_knowledge=10,
                orbit_altitude=820,
                planet='earth')

leo.define_wheel(inertia=100, density=steel_density, thickness_ratio=1/75, inner_outer_ratio=97/100)
print(f'Wheel has radius {leo.wheel_r_outer:.3f} m and weighs {leo.wheel_mass:.3f} kg.')

# gravity gradient disturbance
leo.get_max_h_in_periodic_disturbance()
leo.get_min_h_for_gyro_stiffness()
print(f'Momentum wheel should spin at {ceil(leo.min_omega*60/2/pi)} RPM to attenuate gravity gradient disturbance.')

# manoeuvring requirements
leo.get_slew_torque()