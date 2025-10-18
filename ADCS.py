from math import pi, sqrt, ceil
from main import *

MF = 0.05

class Satellite:
    def __init__(self, inertia, max_ampl_disturbance_torque, thrust_moment_arm, manoeuvre_degrees, manoeuvre_secs, pointing_knowledge, case_name, planet: str):
        self.inertia = inertia # kg m²
        self.max_ampl_disturbance_torque = max_ampl_disturbance_torque # N m
        self.thrust_moment_arm = thrust_moment_arm  # m
        self.manoeuvre_degrees = manoeuvre_degrees # deg
        self.manoeuvre_secs = manoeuvre_secs # s
        self.pointing_knowledge = pointing_knowledge # arcsec
        self.pointing_knowledge_degrees = pointing_knowledge/3600 # deg
        self.orbit_dist = exercise_data[case_name][8] + planet_data[planet]['mean_radius'] # m
        self.orbital_period = 2*pi*sqrt(self.orbit_dist**3/planet_data[planet]['gravitational_parameter'])
        self.max_rpm = 5_000 # https://www.aac-clyde.space/what-we-do/space-products-components/adcs/rw400
        # TODO add MF

    def run(self):
        # manoeuvring requirements
        self.slew_to_rest()
        self.thrust_for_disturbance_counteract()
        self.get_thruster_type(max(self.slew_thrust_to_rest, self.disturbance_counteract_thrust))

        # gravity gradient disturbance
        if self.get_max_undamped_drift() > self.pointing_knowledge_degrees:
            if self.slew_thrust_to_rest > 50*self.disturbance_counteract_thrust and self.thruster_type!='electric':
                print(f'Reaction wheel is needed against disturbance ({self.max_undamped_drift:.3g}>{self.pointing_knowledge_degrees:.3g} deg, and thrust available is {self.slew_thrust_to_rest/self.disturbance_counteract_thrust:.0f}x bigger than needed)')
                self.get_max_h_in_periodic_disturbance_for_reaction_wheel()
                self.define_wheel(self.h_disturb, self.slew_h_to_rest)
                print(f'Reaction wheel will spin at up to {self.h_disturb/self.wheel_inertia*60/2/pi:.0f} RPM counteracting gravity gradient disturbance. ({self.slew_h_to_rest/self.wheel_inertia*60/2/pi:.0f} RPM if used to manoeuvre)')
            else:
                print(f'Disturbance can be counteracted with thrusters: ({self.slew_thrust_to_rest}N vs {self.disturbance_counteract_thrust}N {' but electric propulsion can be throttled low' if self.thruster_type=='electric' else ''})')
        else:
            print(f'Reaction wheel is NOT needed ({self.max_undamped_drift:.3g}<{self.pointing_knowledge_degrees:.3g} deg)')

    def get_max_undamped_drift(self):
        self.max_undamped_drift = (self.orbital_period/2)**2 * self.max_ampl_disturbance_torque / (2 * self.inertia * pi**2) * 180/pi
        return self.max_undamped_drift

    def get_min_h_for_gyro_stiffness_for_momentum_wheel(self):
        '''
        Returns the minimum momentum a wheel needs to have to neutralise periodic disturbance torques to within the allowed range
        '''
        raise SystemError(f'This is not part of the task') 
        self.min_h = self.max_ampl_disturbance_torque * self.orbital_period / (self.pointing_knowledge_degrees*pi/180 * 4) # kg m² / s
        self.min_omega = self.min_h / self.wheel_inertia
        return self.min_h
    
    def get_max_h_in_periodic_disturbance_for_reaction_wheel(self):
        '''
        Returns the maximum momentum a reaction wheel will have in counteracting the periodic disturbance torque
        '''
        self.h_disturb = sqrt(2)/2 * self.max_ampl_disturbance_torque * self.orbital_period/4 # kg m² / s
        return self.h_disturb
    
    def define_wheel(self, angular_momentum1, angular_momentum2):
        inertia1 = angular_momentum1 / (self.max_rpm/60*2*pi) # kg m²
        wheel_dia1 = (32*inertia1/(steel_density*pi*thickness_ratio))**(1/5) # m
        wheel_mass1 = 8 * inertia1 / wheel_dia1**2 # kg
    
        inertia2 = angular_momentum2 / (self.max_rpm/60*2*pi) # kg m²
        wheel_dia2 = (32*inertia2/(steel_density*pi*thickness_ratio))**(1/5) # m
        wheel_mass2 = 8 * inertia2 / wheel_dia2**2 # kg

        if wheel_mass2<=wheel_mass1*5:
            print(f'Reaction wheel can reasonably be used for manoeuvring requirement ({wheel_mass2:.3g}kg vs {wheel_mass1:.3g}kg) (not strictly the task)')
            self.wheel_inertia = inertia2
            self.wheel_dia = wheel_dia2
            self.wheel_mass = wheel_mass2
        else:
            print(f'Reaction wheel cannot be used for manoeuvring requirement ({wheel_mass2:.3g}kg vs {wheel_mass1:.3g}kg)')
            self.wheel_inertia = inertia1
            self.wheel_dia = wheel_dia1
            self.wheel_mass = wheel_mass1
        print(f'Wheel has diameter {self.wheel_dia:.3g} m and weighs {self.wheel_mass:.3g} kg.')

    def slew_to_rest(self):
        '''
        Returns the torque needed to satisfy slew requirements, ending in standstill
        '''
        self.slew_torque_to_rest = 4 * self.inertia * self.manoeuvre_degrees / self.manoeuvre_secs**2 # N m
        self.slew_thrust_to_rest = self.slew_torque_to_rest / self.thrust_moment_arm
        self.slew_h_to_rest = self.slew_torque_to_rest * self.manoeuvre_secs
        return self.slew_torque_to_rest

    def thrust_for_wheel_dump(self, h, time):
        '''
        Returns thrust to dump momentum h in 'time' seconds
        '''
        F = h/(time*self.thrust_moment_arm) # N
        return F
    
    def thrust_for_disturbance_counteract(self):
        self.disturbance_counteract_thrust = self.max_ampl_disturbance_torque * (1+MF) / self.thrust_moment_arm
        return self.disturbance_counteract_thrust

    def get_thruster_type(self, max_force):
        if max_force<0.01:
            self.thruster_type = 'electric'
        elif max_force<1:
            self.thruster_type = 'cold gas'
        elif max_force>1:
            self.thruster_type = 'monopropellant'
        print(f'Thrusters should be {self.thruster_type} (slew: {self.slew_thrust_to_rest:.3g}N, disturbance: {self.disturbance_counteract_thrust:.3g}N).')

steel_density = 8050 # kg/m³
thickness_ratio = 1/8

print('THEOS-2')
theos = Satellite(inertia=200, 
                max_ampl_disturbance_torque=1e-2,
                thrust_moment_arm=1,
                manoeuvre_degrees=25,
                manoeuvre_secs=60,
                pointing_knowledge=10,
                case_name='Case 1',
                planet='earth')
theos.run()

print('\nChang e-4')
change = Satellite(inertia=100, 
                max_ampl_disturbance_torque=1e-4,
                thrust_moment_arm=1,
                manoeuvre_degrees=5,
                manoeuvre_secs=10,
                pointing_knowledge=5,
                case_name='Case 2',
                planet='moon')
change.run()

print('\nMarco')
marco = Satellite(inertia=1e-2, 
                max_ampl_disturbance_torque=1e-7,
                thrust_moment_arm=0.05,
                manoeuvre_degrees=45,
                manoeuvre_secs=60,
                pointing_knowledge=2000,
                case_name='Case 3',
                planet='mars')
marco.run()

print('\nOdyssey')
odyssey = Satellite(inertia=100, 
                max_ampl_disturbance_torque=1e-2,
                thrust_moment_arm=1,
                manoeuvre_degrees=30,
                manoeuvre_secs=30,
                pointing_knowledge=5,
                case_name='Case 4',
                planet='mars')
odyssey.run()

print('\nAkatsuki')
akatsuki = Satellite(inertia=5, 
                max_ampl_disturbance_torque=1e-3,
                thrust_moment_arm=0.5,
                manoeuvre_degrees=10,
                manoeuvre_secs=45,
                pointing_knowledge=5,
                case_name='Case 4',
                planet='venus')
akatsuki.run()

