import tkinter as tk
from tkinter import ttk
import sys
import numpy as np
import matplotlib.pyplot as plt

# Goal: use the link equation in decibel form so Eb/No = sum of all gains and losses
# for each configuration calculate the margin
# note: all losses are maximal losses based on the worst expected conditions
# Assumptions: T0=290K, G/T=G-T0((1-F_loss)/F_loss) [dB]

# constants
k_b = 1.381e-23
c = 3e8
eta_ant = 0.55
min_elev = 10 * np.pi / 180  # rad
EbNo_req = 10.5  # [dB], value comes from BPSK modulation for BER of 10e-6, from ADSEE lecture 4
T_0 = 290
T_ref_sc = 325 #K
T_ref_gs = 60 #K
loss_pointing_gs = 0.12 #dB

# Planetary characteristics database (SI units)

planet_data = {
    "mars": {
        "gravitational_parameter": 4.282837e13,  # m^3/s^2
        "mean_radius": 3.3895e6,  # m
        "orbital_period": 5.935e7,  # s (687 Earth days)
        "rotation_period": 8.861e4,  # s (24.6229 hours)
        "max_distance_to_earth": 4.02e11,  # m
        "min_distance_to_earth": 5.46e10,  # m
        "max_distance_to_sun": 2.492e11  # m
    },
    "venus": {
        "gravitational_parameter": 3.24859e14,  # m^3/s^2
        "mean_radius": 6.0518e6,  # m
        "orbital_period": 1.944e7,  # s (224.7 Earth days)
        "rotation_period": -2.004e6,  # s (retrograde, −243.018 days)
        "max_distance_to_earth": 2.58e11,  # m
        "min_distance_to_earth": 4.2e7,  # m
        "max_distance_to_sun": 1.089e11  # m
    },
    "earth": {
        "gravitational_parameter": 3.986004418e14,  # m^3/s^2
        "mean_radius": 6.371e6,  # m
        "orbital_period": 3.15576e7,  # s (365.25 days)
        "rotation_period": 8.6164e4,  # s (23.934 hours)
        "max_distance_to_earth": 0.0,  # m (reference body)
        "min_distance_to_earth": 0.0,  # m
        "max_distance_to_sun": 1.521e11  # m
    },
    "moon": {
        "gravitational_parameter": 4.9048695e12,  # m^3/s^2
        "mean_radius": 1.7374e6,  # m
        "orbital_period": 2.3606e6,  # s (27.3 Earth days)
        "rotation_period": 2.3606e6,  # s (synchronous rotation)
        "max_distance_to_earth": 4.07e8,  # m
        "min_distance_to_earth": 3.63e8,  # m
        "max_distance_to_sun": 1.521e11  # m  # same as Earth
    },
    "mercury": {
        "gravitational_parameter": 2.2032e13,  # m^3/s^2
        "mean_radius": 2.4397e6,  # m
        "orbital_period": 7.6005e6,  # s (88 Earth days)
        "rotation_period": 5.067e6,  # s (58.646 Earth days)
        "max_distance_to_earth": 2.22e11,  # m
        "min_distance_to_earth": 7.7e10,  # m
        "max_distance_to_sun": 6.981e10  # m
    }
}
exercise_data = {
    'Case 1': [150, 400, 0.8, 0.7, 2.2, 0.920833, 1, 15, 820, 0, 0.12, 100, 45, 0.01, 32, 100, 0.5, 10.5, 0, 'earth'],
    'Case 2': [100, 400, 0.8, 0.7, 2.2, 0.920833, 4.2, 15, 100, 0, 0.1, 10, 45, 0.05, 8, 75, 8, 10.5, 0, 'moon'],
    'Case 3': [10, 400, 0.8, 0.7, 8.4, 0.851136, 0.1, 35, 1000, 5, 1, 0.01, 10, 0.2, 8, 5, 18, 10.5,1, 'mars'],
    'Case 4': [100, 1000, 0.8, 0.7, 8.4, 0.851136, 2, 35, 400, 20, 0.1, 1, 10, 0.1, 8, 10, 18, 10.5,1, 'mars'],
    'Case 5': [200, 1000, 0.8, 0.7, 8.5, 0.851136, 1.6, 35, 1000, 10, 0.1, 1, 10, 0.3, 8, 25, 12, 10.5,1, 'venus']
}

# formulas
def loss_space(values, f, c, planet_data, min_elev, typ):
    # procedure taken from SMAD 3rd edition page 113
    if values[-1] == "earth":
        Re = planet_data["earth"]['mean_radius']
        h = values[8] * 1e3
        sin_rho = Re / (Re + h)
        eta = np.arcsin(sin_rho * np.cos(min_elev))
        l = (np.pi / 2) - min_elev - eta
        d = Re * (np.sin(l) / np.sin(eta))
    elif values[-1] == "moon":
        Re = planet_data["earth"]['mean_radius']
        h = planet_data["moon"]['max_distance_to_earth']
        sin_rho = Re / (Re + h)
        eta = np.arcsin(sin_rho * np.cos(min_elev))
        l = np.pi / 2 - min_elev - eta
        d = Re * (np.sin(l) / np.sin(eta))
    else:
        if typ == 1:
            d_e = planet_data[values[-1]]['max_distance_to_earth']
            d_s = planet_data[values[-1]]['max_distance_to_sun']
            d = np.sqrt(d_e ** 2 + d_s ** 2 - 2 * d_e * d_s * np.cos(values[9]*np.pi/180))
        else:
            d = planet_data[values[-1]]['max_distance_to_earth']
    loss_space = 20 * np.log10((4 * np.pi * d * f) / c)
    return -loss_space
def eirp(values, f, D, P, c, eta_ant):
    loss_tx = 10 * np.log10(values[2])
    gain = 10 * np.log10(eta_ant * ((np.pi * D) / (c / f)) ** 2)
    eirp = 10 * np.log10(P) + gain + loss_tx
    return eirp
def g_over_t(values, f, D, T_0,T_ref, c, eta_ant):
    T_sys = 10 * np.log10(T_ref+T_0 * (1 - values[3])/values[3])
    gain = 10 * np.log10(eta_ant * ((np.pi * D) / (c / f))**2)
    g_over_t = gain -T_sys
    return g_over_t
def loss_atm(f, min_elev):
    # linear approximation of SMAD 3rd edition book p565
    f_GHz = f/1e9
    loss_atm = 0.04 + 0.002 * (f_GHz / np.sin(min_elev))
    return -loss_atm
def loss_pointing(f,D,pointing_accuracy,mode):
    f_GHz = f/1e9
    half_power_beamwidth = 21 / (f_GHz * D)
    if mode == 1:
        pointing_accuracy=0.1*half_power_beamwidth
    loss_pointing = 12 * (pointing_accuracy / half_power_beamwidth) ** 2
    return -loss_pointing
def bitrate_loss(datarate):
    bitrate_loss = 10 * np.log10(datarate)
    return -bitrate_loss
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
    f_downlink = values[4]*1e9
    f_uplink = values[4] * values[5] * 1e9
    uplink_datarate = values[11] * 1e6
    downlink_datarate_val = downlinkdatarate(values, planet_data)
    typ = values[18]

    EbNo_downlink_values = [
            eirp(values, f_downlink, values[6], values[0], c, eta_ant),
            g_over_t(values, f_downlink, values[7], T_0, T_ref_gs, c, eta_ant),
            bitrate_loss(downlink_datarate_val),
            -10*np.log10(k_b),
            loss_space(values, f_downlink, c, planet_data, min_elev,typ),
            loss_atm(f_downlink, min_elev),
            loss_pointing(f_downlink, values[6],values[10],0)
    ]
    EbNo_uplink_values = [
            eirp(values, f_uplink, values[7], values[1], c, eta_ant),
            g_over_t(values, f_uplink, values[6], T_0, T_ref_sc, c, eta_ant),
            bitrate_loss(uplink_datarate),
            -10 * np.log10(k_b),
            loss_space(values, f_uplink, c, planet_data, min_elev,typ),
            loss_atm(f_uplink, min_elev),
            loss_pointing(f_uplink, values[7],0,1)
    ]
    margin_downlink = sum(EbNo_downlink_values) - values[17]
    margin_uplink = sum(EbNo_uplink_values) - values[17]

    return margin_downlink, margin_uplink, EbNo_downlink_values, EbNo_uplink_values

def main():

    # GUI
    def submit():
        raw_values = [entry.get() for entry in entries]
        values = []
        for value in raw_values[:-1]:
            if not value:
                value = 0.0
            values.append(float(value))
        values.append(raw_values[-1].lower())

        margin_downlink, margin_uplink, downlink_values, uplink_values = link(values, eta_ant, c, k_b, min_elev)

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
    def choose_mode():
        # Create a top-level window
        window = tk.Tk()
        window.title("Select Mode")
        window.geometry("300x100")
        window.resizable(False, False)

        # Logic variable to store result
        selected_mode = tk.StringVar(value="")  # default empty

        def set_mode(mode):
            selected_mode.set(mode)
            window.quit()
            window.destroy()

        # UI elements
        tk.Label(window, text="Choose mode:", font=("Arial", 12)).pack(pady=10)

        button_frame = tk.Frame(window)
        button_frame.pack()

        tk.Button(button_frame, text="Custom", width=10, command=lambda: set_mode("custom")).pack(side="left", padx=10)
        tk.Button(button_frame, text="Exercise", width=10, command=lambda: set_mode("exercise")).pack(side="right", padx=10)

        # Start GUI loop
        window.mainloop()

        return selected_mode.get()
    def show_margin_popup(margin_cases):
        popup = tk.Toplevel()
        popup.title("Link Margin Results")
        popup.geometry("1300x450")
        popup.grab_set()

        columns = ("Case", "Link Type","Required Eb/No", "Margin", "EIRP", "G/T", "Bitrate",
                "Free Space Loss", "Atmospheric Loss", "Pointing Loss")

        tree = ttk.Treeview(popup, columns=columns, show="headings", height=18)

        for col in columns:
            tree.heading(col, text=col)
            tree.column(col, anchor="center", width=120)

        for i in range(1, 6):
            # Downlink row
            downlink = f"{margin_cases[i][0]:.2f} dB"
            down_data = margin_cases[i][2]
            tree.insert("", "end", values=(
                f"Case {i}", "Downlink", f"{exercise_data[f'Case {i}'][17]:.2f}", downlink,
                f"{down_data[0]:.2f} dB" , f"{down_data[1]:.2f} dB", f"{down_data[2]:.2f} dB",
                f"{down_data[4]:.2f} dB", f"{down_data[5]:.2f} dB", f"{down_data[6]:.2f} dB"
            ))

            # Uplink row
            uplink = f"{margin_cases[i][1]:.2f} dB"
            up_data = margin_cases[i][3]
            tree.insert("", "end", values=(
                f"Case {i}", "Uplink", f"{exercise_data[f'Case {i}'][17]:.2f}",uplink,
                f"{up_data[0]:.2f} dB", f"{up_data[1]:.2f} dB", f"{up_data[2]:.2f} dB",
                f"{up_data[4]:.2f} dB", f"{up_data[5]:.2f} dB", f"{up_data[6]:.2f} dB"
            ))

            # Divider row
            tree.insert("", "end", values=("────────────", "────────────", "────────────", "────────────", "────────────", "────────────", "────────────", "────────────", "────────────", "────────────"))

        tree.pack(fill="both", expand=True, padx=10, pady=10)

        def on_ok():
            popup.destroy()
            sys.exit()

        ok_button = ttk.Button(popup, text="OK", command=on_ok)
        ok_button.pack(pady=5)

    selected_mode = choose_mode()
    if selected_mode == "custom":
        root = tk.Tk()
        root.title("Enter Transmission Parameters")
        root.geometry("500x800")
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
            "13. Payload pixel size [arcmin]",
            "14. Payload bits per pixel [bit]",
            "15. Payload duty cycle [%]",
            "16. Payload downlink time [hr]",
            "17. Required Eb/No [dB]",
            "18. Type(default max distance, 1 for elongation angle)",
            "19. Planet",
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
    elif selected_mode == "exercise":
        margin_cases = {
            i: link(exercise_data[f'Case {i}'], eta_ant, c, k_b, min_elev)
            for i in range(1, 6)
        }
        root = tk.Tk()
        root.withdraw()
        show_margin_popup(margin_cases)
        root.mainloop()

if __name__=='__main__':
    main()