import tkinter as tk
from traceback import print_tb

import numpy as np
import matplotlib.pyplot as plt

# Goal: use the link equation in decibel form so Eb/No = sum of all gains and losses
# for each configuration calculate the margin

#constants
k_b = 1.381*10e-23
c = 10**8

#formulas
def Loss_space(values,d,f,c):
    d = values[8]
    f = values[9]
    loss_space = 20*np.log10((4*np.pi*d)/(c/f))
    return loss_space


#GUI
def submit():
    values = [float(entry.get()) for entry in entries]
    print(values)

root = tk.Tk()
root.title("Enter Transmission Parameters")
root.geometry("400x600")

labels = [
    "Transmitter power (spacecraft) [W]",
    "Transmitter power (ground station) [W]",
    "Loss factor transmitter [dB]",
    "Loss factor receiver [dB]",
    "Downlink frequency [Hz]",
    "Turnaround ratio (uplink/downlink) [-]",
    "Antenna diameter (spacecraft) [m]",
    "Antenna diameter (ground station) [m]",
    "Orbit altitude [km]",
    "GS Elevation angle [deg]",
    "Pointing offset angle [deg]",
    "Required uplink data rate [bit/s]",
    "Payload swath width angle [deg]",
    "Payload pixel size [m]",
    "Payload bits per pixel [bit]",
    "Payload duty cycle [%]",
    "Required data rate [bit/s]"
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
