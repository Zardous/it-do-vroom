import tkinter as tk
from traceback import print_tb

import numpy as np
import matplotlib.pyplot as plt

# Goal: use the link equation in decibel form so Eb/No = sum of all gains and losses
# for each configuration calculate the margin

#constants
k_b = 1.381*10e-23
c = 3*10**8
eta_ant = 0.55

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
    }
}

#formulas
def loss_space(values,c,planet_data):
    if values[-1] == "earth":
        l = np.arccos(planet_data["earth"]["mean_radius"]/(planet_data["earth"]["mean_radius"]+values[8]))
        d = np.sin(l)*(planet_data["earth"]["mean_radius"]+values[8])
    else:
        d = planet_data[values[-1]]["min_distance_to_earth"]
    f = values[4]*10**9
    loss_space = 20*np.log10((4*np.pi*d)/(c/f))
    return loss_space

def gain_sc(values, eta_ant,c):
    f = values[4]*10**9
    D = values[6]
    gain_sc = 10*np.log10(eta_ant*((np.pi*D)/(c/f))**2)
    return gain_sc

def gain_gs(values, eta_ant,c):
    f = values[4]*10**9
    D = values[7]
    gain_gs = 10 * np.log10(eta_ant*((np.pi*D)/(c/f))**2)
    return gain_gs

def loss_tx(values):
    loss_tx = 10*np.log10(values[2])
    return loss_tx

def loss_rx(values):
    loss_rx = 10*np.log10(values[3])
    return loss_rx

def eirp_sc(values, gain, loss_tx):
    eirp_sc = 10*np.log10(values[0])+ gain_sc(values,eta_ant,c) - loss_tx(values)
    return eirp_sc

def eirp_gs(values, gain, loss_tx):
    eirp_gs = 10*np.log10(values[1])+ gain_gs(values,eta_ant,c) - loss_tx(values)
    return eirp_gs

#GUI
def submit():
    raw_values = [entry.get() for entry in entries]
    values = []
    for value in raw_values[:-1]:
        if not value:
            value = 0.0
        values.append(float(value))
    values.append(raw_values[-1].lower())
    print(values)
    print(loss_space(values,c,planet_data))
    return values

root = tk.Tk()
root.title("Enter Transmission Parameters")
root.geometry("400x700")

labels = [
    "Transmitter power (spacecraft) [W]",
    "Transmitter power (ground station) [W]",
    "Loss factor transmitter [-]",
    "Loss factor receiver [-]",
    "Downlink frequency [GHz]",
    "Turnaround ratio (uplink/downlink) [-]",
    "Antenna diameter (spacecraft) [m]",
    "Antenna diameter (ground station) [m]",
    "Orbit altitude [km]",
    "Elongation angle [deg]",
    "Pointing offset angle [deg]",
    "Required uplink data rate [bit/s]",
    "Payload swath width angle [deg]",
    "Payload pixel size [m]",
    "Payload bits per pixel [bit]",
    "Payload duty cycle [%]",
    "Required data rate [bit/s]",
    "Planet"
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

