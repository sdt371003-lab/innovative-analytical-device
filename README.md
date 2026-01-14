import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --------------------------------------------------
# Page setup
# --------------------------------------------------
st.set_page_config(page_title="GC-MS vs Integrated System Simulation", layout="wide")
st.title("Simulation: Integrated Analytical System vs Traditional GC-MS")

# --------------------------------------------------
# System selection
# --------------------------------------------------
system_choice = st.radio(
    "Select system to simulate:",
    ("Integrated System (Your Device)", "Traditional GC-MS")
)

# --------------------------------------------------
# Sample input
# --------------------------------------------------
st.header("Sample Input")

compound_names = st.text_area(
    "Enter compound names (comma separated):",
    "Acetone, Methanol, Ethanol"
)
boiling_points = st.text_area(
    "Enter boiling points °C (same order):",
    "56, 65, 78"
)
ion_eff = st.text_area(
    "Enter ionization efficiency % (same order):",
    "70, 60, 65"
)
polarity_factor = st.text_area(
    "Enter polarity factor (0–1) for each compound (optional):",
    "0.2, 0.3, 0.5"
)

# Clean inputs
compound_names = [c.strip() for c in compound_names.split(",") if c.strip()]
boiling_points = [float(b.strip()) for b in boiling_points.split(",") if b.strip()]
ion_eff = [float(i.strip()) / 100 for i in ion_eff.split(",") if i.strip()]
polarity_factor = [float(p.strip()) for p in polarity_factor.split(",") if p.strip()]

if not (len(compound_names) == len(boiling_points) == len(ion_eff) == len(polarity_factor)):
    st.error("⚠️ Please enter the same number of values for all fields.")
    st.stop()

n_compounds = len(compound_names)

# --------------------------------------------------
# Thermal & system settings
# --------------------------------------------------
st.header("System Settings")

if system_choice == "Integrated System (Your Device)":
    # Sample Prep Oven
    prep_temp = st.number_input("Sample Prep Oven Temperature (°C)", 50, 400, 150)
    prep_pressure = st.number_input("Sample Prep Oven Pressure (kPa)", 10.0, 500.0, 101.3)

    # µTD + Hot Channel
    hot_channel_temp = st.number_input("Hot Channel / µTD Temperature (°C)", 100, 400, 250)

    # Column properties
    column_rate = st.number_input("Column heating rate (°C/min)", 1, 50, 5)
    column_length = st.number_input("Column length factor (arbitrary, affects RT & width)", 1.0, 5.0, 1.0)
    stationary_phase_factor = st.number_input("Stationary phase factor (0.5–1.5)", 0.5, 1.5, 1.0)
    detector_sensitivity = st.number_input("Detector sensitivity factor (0.5–2.0)", 0.5, 2.0, 1.0)

    # Hot-factor: focusing & safe transfer
    hot_factor = 1 + (hot_channel_temp - 150)/400 + (column_rate/50)*0.1

else:
    # Traditional GC-MS
    column_rate = st.number_input("Column heating rate (°C/min)", 1, 50, 5)
    column_length = st.number_input("Column length factor (arbitrary, affects RT & width)", 1.0, 5.0, 1.0)
    stationary_phase_factor = st.number_input("Stationary phase factor (0.5–1.5)", 0.5, 1.5, 1.0)
    detector_sensitivity = st.number_input("Detector sensitivity factor (0.5–2.0)", 0.5, 2.0, 1.0)

# --------------------------------------------------
# Time axis
# --------------------------------------------------
x = np.linspace(0, 20, 1200)

# --------------------------------------------------
# Signal containers
# --------------------------------------------------
y_column_ion = np.zeros_like(x)
y_column_neutral = np.zeros_like(x)
peak_table = []

# --------------------------------------------------
# Resolution function
# --------------------------------------------------
def resolution(t1, w1, t2, w2):
    return 2 * abs(t2 - t1) / (w1 + w2)

# --------------------------------------------------
# Simulation loop
# --------------------------------------------------
for i in range(n_compounds):
    # Retention time (approx.): boiling point + polarity + column factors
    center = 2 + i*2 + (boiling_points[i]/100)*stationary_phase_factor*column_length + polarity_factor[i]*1.0

    # Peak width (approx.)
    width = 0.15 + (boiling_points[i]/500)*(1/column_rate)*column_length

    if system_choice == "Integrated System (Your Device)":
        # Hot channel / µTD improves focusing
        width_column = width / hot_factor
        intensity_ion = ion_eff[i]*hot_factor*detector_sensitivity
        intensity_neutral = (1-ion_eff[i])*hot_factor*detector_sensitivity
    else:
        width_column = width
        intensity_ion = ion_eff[i]*detector_sensitivity
        intensity_neutral = 0

    # Gaussian peaks
    peak_ion = intensity_ion * np.exp(-((x - center)**2) / (2*width_column**2))
    peak_neutral = intensity_neutral * np.exp(-((x - center)**2) / (2*width_column**2))

    y_column_ion += peak_ion
    y_column_neutral += peak_neutral

    peak_table.append({
        "Compound": compound_names[i],
        "Retention Time (min)": round(center,2),
        "Peak Width": round(width_column,3),
        "Ionic Intensity": round(intensity_ion,3),
        "Neutral Intensity": round(intensity_neutral,3)
    })

# --------------------------------------------------
# Resolution calculation
# --------------------------------------------------
for i in range(1, len(peak_table)):
    t1 = peak_table[i-1]["Retention Time (min)"]
    t2 = peak_table[i]["Retention Time (min)"]
    w1 = peak_table[i-1]["Peak Width"]
    w2 = peak_table[i]["Peak Width"]
    peak_table[i]["Resolution vs Previous"] = round(resolution(t1,w1,t2,w2),2)
peak_table[0]["Resolution vs Previous"] = None

df = pd.DataFrame(peak_table)

# --------------------------------------------------
# Plot
# --------------------------------------------------
st.header("Simulated Chromatogram")
fig, ax = plt.subplots(figsize=(14,6))
ax.plot(x, y_column_ion, color="blue", label="Column – Ionic Path")
if system_choice == "Integrated System (Your Device)":
    ax.plot(x, y_column_neutral, color="green", linestyle="--", label="Column – Neutral Path")

for row in peak_table:
    ax.scatter(row["Retention Time (min)"], row["Ionic Intensity"], zorder=3, color="blue")
    ax.text(row["Retention Time (min)"], row["Ionic Intensity"]+0.02, row["Compound"], ha="center")

ax.set_xlabel("Retention Time (min)")
ax.set_ylabel("Signal (a.u.)")
ax.set_title(system_choice)
ax.legend()
st.pyplot(fig)

# --------------------------------------------------
# Results table
# --------------------------------------------------
st.header("Detected Peaks & Resolution")
st.dataframe(df)

# --------------------------------------------------
# Option to download table
# --------------------------------------------------
csv = df.to_csv(index=False).encode()
st.download_button(
    label="Download peaks table as CSV",
    data=csv,
    file_name='simulated_peaks.csv',
    mime='text/csv'
)

# --------------------------------------------------
# Scientific notes
# --------------------------------------------------
st.info(
    "🔹 This simulation is educational.\n"
    "Sample Prep Oven controls sample conditioning (temperature & pressure) and does NOT directly affect peak shapes.\n"
    "Hot Channel & µTD provide safe transfer and sample focusing, improving separation and sample integrity.\n"
    "Retention times, peak widths, and intensities are approximations; real GC-MS behavior depends on column chemistry, length, carrier gas, pressure, and detector physics."
)
