import tkinter as tk
from traceback import print_tb

import numpy as np
import matplotlib.pyplot as plt

# Goal: use the link equation in decibel form so Eb/No = sum of all gains and losses
# for each configuration calculate the margin
# note: all losses are maximal losses based on the worst expected conditions
#Assumptions: T0=290K, G/T=G-T0(1-F_loss) [dB]

#constants
k_b = 1.381*10e-23
c = 3*10**8
eta_ant = 0.55
min_elev = 10* np.pi / 180 #rad
EbNo_req = 10.5 #[dB], value comes from BPSK modulation for BER of 10e-6, from ADSEE lecture 4

# Planetary characteristics database (SI units)

planet_data = {
    "mars": {
        "gravitational_parameter": 4.282837e13,       # m^3/s^2
        "mean_radius": 3.3895e6,                      # m
        "orbital_period": 5.935e7,                    # s (687 Earth days)
        "rotation_period": 8.861e4,                   # s (24.6229 hours)
        "max_distance_to_earth": 4.02e11,             # m
        "min_distance_to_earth": 5.46e10              # m
    },
    "venus": {
        "gravitational_parameter": 3.24859e14,        # m^3/s^2
        "mean_radius": 6.0518e6,                      # m
        "orbital_period": 1.944e7,                    # s (224.7 Earth days)
        "rotation_period": -2.004e6,                  # s (retrograde, −243.018 days)
        "max_distance_to_earth": 2.58e11,             # m
        "min_distance_to_earth": 4.2e7                # m
    },
    "earth": {
        "gravitational_parameter": 3.986004418e14,    # m^3/s^2
        "mean_radius": 6.371e6,                       # m
        "orbital_period": 3.15576e7,                  # s (365.25 days)
        "rotation_period": 8.6164e4,                  # s (23.934 hours)
        "max_distance_to_earth": 0.0,                 # m (reference body)
        "min_distance_to_earth": 0.0                  # m
    },
    "moon": {
        "gravitational_parameter": 4.9048695e12,      # m^3/s^2
        "mean_radius": 1.7374e6,                      # m
        "orbital_period": 2.3606e6,                   # s (27.3 Earth days)
        "rotation_period": 2.3606e6,                  # s (synchronous rotation)
        "max_distance_to_earth": 4.07e8,              # m
        "min_distance_to_earth": 3.63e8               # m
    },
    "mercury": {
            "gravitational_parameter": 2.2032e13,         # m^3/s^2
            "mean_radius": 2.4397e6,                      # m
            "orbital_period": 7.6005e6,                   # s (88 Earth days)
            "rotation_period": 5.067e6,                   # s (58.646 Earth days)
            "max_distance_to_earth": 2.22e11,             # m
            "min_distance_to_earth": 7.7e10               # m
    }
}

#formulas
def loss_space(values,f,c,planet_data,min_elev):
    #procedure taken from SMAD 3rd edition page 113
    if values[-1] == "earth":
        Re = planet_data["earth"]['mean_radius']
        h = values[8]*1000
        sin_rho = Re/(Re+h)
        eta = np.arcsin(sin_rho*np.cos(min_elev))
        l = np.pi/2 - min_elev - eta
        d = Re*(np.sin(l)/np.sin(eta))
    else:
        d = planet_data[values[-1]]["max_distance_to_earth"]
    loss_space = 20*np.log10((4*np.pi*d)/(c/f))
    return loss_space
def eirp(values,f,D,P,c,eta_ant):
    loss_tx = 10 * np.log10(values[2])
    gain = 10 * np.log10(eta_ant * ((np.pi * D) / (c / f)) ** 2)
    eirp = 10*np.log10(P)+ gain - loss_tx
    return eirp
def g_over_t(values,f,D,c,eta_ant):
    T_sys = 10 * np.log10(290*(1-values[3]))
    gain = 10 * np.log10(eta_ant * ((np.pi * D) / (c / f)) ** 2)
    g_over_t = gain - T_sys
    return g_over_t
def loss_atm(f,min_elev):
    #linear approximation of SMAD 3rd edition book p565
    f_GHz = f/(10**9)
    loss_atm = 0.04+0.002*f_GHz/np.sin(min_elev)
    return loss_atm
def loss_pointing(values):
    beamwidth = 21/(values[4]*values[6])
    loss_pointing = 12*(values[10]/beamwidth)**2
    return loss_pointing
def bitrate_term_loss(datarate,k_b):
    bitrate_term_loss = 10*np.log10(datarate*k_b)
    return bitrate_term_loss
def downlinkdatarate(values, planet_data):
    rad_planet = planet_data[values[-1]]['mean_radius']
    h_orbit = values[8]
    gravparam = planet_data[values[-1]]['gravitational_parameter']
    swathwidth = values[12]
    pixelsize = values[13]
    bitsperpixel = values[14]
    v_orb = np.sqrt(gravparam / (rad_planet + h_orbit))
    v_ground = v_orb * rad_planet / (rad_planet + h_orbit)
    swath_time = h_orbit * np.tan(pixelsize) / v_ground
    pixelsperswath = swathwidth / pixelsize
    pixelrate = pixelsperswath / swath_time
    payloaddatarate = pixelrate * bitsperpixel
    downlinktime = values[17]
    dutycycle = values[15]
    downlinkdatarate = 24 * 60 * 60 * dutycycle / 100 * payloaddatarate / (downlinktime * 60 * 60)
    return downlinkdatarate
def link(values, eta_ant, c, k_b, min_elev):

    f_downlink = values[4] * 10**9
    f_uplink = values[4] * values[5] * 10**9
    uplink_datarate = values[11]*10**6
    downlink_datarate_val = downlinkdatarate(values, planet_data)

    EbNo_downlink = (
        eirp(values, f_downlink, values[6],values[0], c, eta_ant)
        + g_over_t(values, f_downlink, values[7], c, eta_ant)
        - bitrate_term_loss(downlink_datarate_val, k_b)
        - loss_space(values, f_downlink, c, planet_data, min_elev)
        - loss_atm(f_downlink, min_elev)
        - loss_pointing(values)
    )
    EbNo_uplink = (
        eirp(values, f_uplink, values[7], values[1] , c, eta_ant)
        + g_over_t(values, f_uplink, values[6], c, eta_ant)
        - bitrate_term_loss(uplink_datarate, k_b)
        - loss_space(values, f_uplink, c, planet_data, min_elev)
        - loss_atm(f_uplink, min_elev)
        - loss_pointing(values)
    )

    return EbNo_downlink, EbNo_uplink

#GUI
def submit():
    raw_values = [entry.get() for entry in entries]
    values = []
    for value in raw_values[:-1]:
        if not value:
            value = 0.0
        values.append(float(value))
    values.append(raw_values[-1].lower())

    EbNo_downlink, EbNo_uplink = link(values, eta_ant, c, k_b, min_elev)
    EbNo_required = values[17]
    margin_downlink, margin_uplink = EbNo_downlink - EbNo_required, EbNo_uplink - EbNo_required

    # Set color based on margin value
    if margin_downlink < 0:
        color_down = "red"
    elif margin_downlink < 3:
        color_down = "orange"
    else:
        color_down = "green"

    if margin_uplink < 0:
        color_up = "red"
    elif margin_uplink < 3:
        color_up = "orange"
    else:
        color_up = "green"

    margin_label_down.config(
        text=f"Downlink Margin: {margin_downlink:.2f} dB",
        fg=color_down
    )

    margin_label_up.config(
        text=f"Uplink Margin: {margin_uplink:.2f} dB",
        fg=color_up
    )

    return values


root = tk.Tk()
root.title("Enter Transmission Parameters")
root.geometry("500x750")
labels = [
    "0. Transmitter power (spacecraft) [W]",
    "1. Transmitter power (ground station) [W]",
    "2. Loss factor transmitter [-]",
    "3. Loss factor receiver [-]",
    "4. Downlink frequency [GHz]",
    "5. Turnaround ratio (uplink/downlink) [-]",
    "6. Antenna diameter (spacecraft) [m]",
    "7. Antenna diameter (ground station) [m]",
    "8. Orbit altitude [km]",
    "9. Elongation angle [deg]",
    "10. Pointing offset angle [deg]",
    "11. Required uplink data rate [Mbit/s]",
    "12. Payload swath width angle [deg]",
    "13. Payload pixel size [m]",
    "14. Payload bits per pixel [bit]",
    "15. Payload duty cycle [%]",
    "16. Payload downlink time [hr]",
    "17. Required Eb/No [dB]",
    "18. Planet"
]
entries = []
for i, label_text in enumerate(labels):
    label = tk.Label(root, text=label_text, anchor="w")
    label.grid(row=i, column=0, padx=10, pady=5, sticky="w")
    entry = tk.Entry(root, width=15)
    entry.grid(row=i, column=1, padx=10, pady=5)
    entries.append(entry)

submit_btn = tk.Button(root, text="Submit", command=submit)
submit_btn.grid(row=len(labels), column=0, columnspan=2, pady=15)
margin_label_down = tk.Label(root, text="", font=("Arial", 12, "bold"))
margin_label_down.grid(row=len(labels) + 1, column=0, columnspan=2, pady=5)

margin_label_up = tk.Label(root, text="", font=("Arial", 12, "bold"))
margin_label_up.grid(row=len(labels) + 2, column=0, columnspan=2, pady=5)


root.mainloop()

