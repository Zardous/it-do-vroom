import tkinter as tk
from traceback import print_tb

import numpy as np
import matplotlib.pyplot as plt

from tempcharlotte import downlinkdatarate

# Goal: use the link equation in decibel form so Eb/No = sum of all gains and losses
# for each configuration calculate the margin
# note: all losses are maximal losses based on the worst expected conditions

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
        d = planet_data[values[-1]]["min_distance_to_earth"]
    loss_space = 20*np.log10((4*np.pi*d)/(c/f))
    return loss_space
def eirp(values,f,D,c,eta_ant):
    loss_tx = 10 * np.log10(values[2])
    gain = 10 * np.log10(eta_ant * ((np.pi * D) / (c / f)) ** 2)
    eirp = 10*np.log10(values[0])+ gain - loss_tx
    return eirp
def g_over_t(values,f,D,c,eta_ant):
    loss_rx = 10 * np.log10(values[3])
    gain = 10 * np.log10(eta_ant * ((np.pi * D) / (c / f)) ** 2)
    g_over_t = gain - loss_rx
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
    datarate = 10**8 #for testing purposes
    bitrate_term_loss = 10*np.log10(datarate*k_b)
    return bitrate_term_loss
def link(values, eta_ant, c, k_b, min_elev):
    f_downlink = values[4] * 10**9
    f_uplink = values[4] * values[5] * 10**9

    # Calculate each term individually
    eirp_val = eirp(values, f_downlink, values[6], c, eta_ant)
    g_over_t_val = g_over_t(values, f_downlink, values[7], c, eta_ant)
    bitrate_loss_val = bitrate_term_loss(values[16], k_b)
    space_loss_val = loss_space(values, f_downlink, c, planet_data, min_elev)
    atm_loss_val = loss_atm(f_downlink, min_elev)
    pointing_loss_val = loss_pointing(values)

    # DEBUG
    print(f"EIRP: {eirp_val:.2f} dB")
    print(f"G/T: {g_over_t_val:.2f} dB/K")
    print(f"Bitrate Term Loss: {bitrate_loss_val:.2f} dB")
    print(f"Free-Space Path Loss: {space_loss_val:.2f} dB")
    print(f"Atmospheric Loss: {atm_loss_val:.2f} dB")
    print(f"Pointing Loss: {pointing_loss_val:.2f} dB")

    EbNo_downlink = (
        eirp_val
        + g_over_t_val
        - bitrate_loss_val
        - space_loss_val
        - atm_loss_val
        - pointing_loss_val
    )
    return EbNo_downlink

#GUI
def submit():
    raw_values = [entry.get() for entry in entries]
    values = []
    for value in raw_values[:-1]:
        if not value:
            value = 0.0
        values.append(float(value))
    values.append(raw_values[-1].lower())
    print(f"Eb/No (Downlink): {link(values,eta_ant,c,k_b,min_elev):.2f} dB")
    print(f"Margin :{link(values,eta_ant,c,k_b,min_elev)-values[16]:.2f} dB")
    return values

root = tk.Tk()
root.title("Enter Transmission Parameters")
root.geometry("400x700")

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
    "11. Required uplink data rate [bit/s]",
    "12. Payload swath width angle [deg]",
    "13. Payload pixel size [m]",
    "14. Payload bits per pixel [bit]",
    "15. Payload duty cycle [%]",
    "16. Required Eb/no [-]",
    "17. Planet"
]

entries = []
for i, label_text in enumerate(labels):
    label = tk.Label(root, text=label_text, anchor="w")
    label.grid(row=i, column=0, padx=10, pady=5, sticky="w")
    entry = tk.Entry(root, width=30)
    entry.grid(row=i, column=1, padx=10, pady=5)
    entries.append(entry)

submit_btn = tk.Button(root, text="Submit", command=submit)
submit_btn.grid(row=len(labels), column=0, columnspan=2, pady=15)

root.mainloop()

