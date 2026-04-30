import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

WATER = "Water"

# -----------------------------
# Property functions
# -----------------------------
def h_PT(P, T):
    return CP.PropsSI("H", "P", P * 1000, "T", T + 273.15, WATER) / 1000

def s_PT(P, T):
    return CP.PropsSI("S", "P", P * 1000, "T", T + 273.15, WATER) / 1000

def h_Ps(P, s):
    return CP.PropsSI("H", "P", P * 1000, "S", s * 1000, WATER) / 1000

def hf(P):
    return CP.PropsSI("H", "P", P * 1000, "Q", 0, WATER) / 1000

def sf(P):
    return CP.PropsSI("S", "P", P * 1000, "Q", 0, WATER) / 1000

def vf(P):
    rho = CP.PropsSI("D", "P", P * 1000, "Q", 0, WATER)
    return 1 / rho

def T_sat(P):
    return CP.PropsSI("T", "P", P * 1000, "Q", 0, WATER) - 273.15

def T_Ps(P, s):
    return CP.PropsSI("T", "P", P * 1000, "S", s * 1000, WATER) - 273.15

def pump(h, v, P1, P2):
    return h + v * (P2 - P1)


# -----------------------------
# Main cycle solver
# -----------------------------
def solve(Pb=8000, Pc=10, Pr=300, Ph=6000, Pl=3500, Po=100, Tin=400, m_dot_boiler=1.0):
    stt = {}

    # Feedwater path
    stt[1] = {"P_kPa": Pc, "T_C": T_sat(Pc), "h": hf(Pc), "s": sf(Pc), "v": vf(Pc)}
    stt[2] = {"P_kPa": Po, "T_C": T_sat(Po), "h": pump(stt[1]["h"], stt[1]["v"], Pc, Po), "s": None, "v": stt[1]["v"]}

    stt[3] = {"P_kPa": Po, "T_C": T_sat(Po), "h": hf(Po), "s": sf(Po), "v": vf(Po)}
    stt[4] = {"P_kPa": Pl, "T_C": T_sat(Pl), "h": pump(stt[3]["h"], stt[3]["v"], Po, Pl), "s": None, "v": stt[3]["v"]}

    stt[5] = {"P_kPa": Pl, "T_C": T_sat(Pl), "h": hf(Pl), "s": sf(Pl), "v": vf(Pl)}
    stt[6] = {"P_kPa": Ph, "T_C": T_sat(Ph), "h": pump(stt[5]["h"], stt[5]["v"], Pl, Ph), "s": None, "v": stt[5]["v"]}

    # Low closed FWH drain / stream to mixer
    stt[7] = {"P_kPa": Pl, "T_C": T_sat(Pl), "h": hf(Pl), "s": sf(Pl), "v": vf(Pl)}

    # Turbine states
    stt[13] = {"P_kPa": Pb, "T_C": Tin, "h": h_PT(Pb, Tin), "s": s_PT(Pb, Tin), "v": None}
    stt[14] = {"P_kPa": Ph, "T_C": T_Ps(Ph, stt[13]["s"]), "h": h_Ps(Ph, stt[13]["s"]), "s": stt[13]["s"], "v": None}
    stt[15] = {"P_kPa": Pl, "T_C": T_Ps(Pl, stt[13]["s"]), "h": h_Ps(Pl, stt[13]["s"]), "s": stt[13]["s"], "v": None}
    stt[16] = {"P_kPa": Pr, "T_C": T_Ps(Pr, stt[13]["s"]), "h": h_Ps(Pr, stt[13]["s"]), "s": stt[13]["s"], "v": None}

    stt[17] = {"P_kPa": Pr, "T_C": Tin, "h": h_PT(Pr, Tin), "s": s_PT(Pr, Tin), "v": None}
    stt[18] = {"P_kPa": Po, "T_C": T_Ps(Po, stt[17]["s"]), "h": h_Ps(Po, stt[17]["s"]), "s": stt[17]["s"], "v": None}
    stt[19] = {"P_kPa": Pc, "T_C": T_Ps(Pc, stt[17]["s"]), "h": h_Ps(Pc, stt[17]["s"]), "s": stt[17]["s"], "v": None}

    # Iterative bleed fractions because state 8 depends on yL and yH
    yH = 0.05
    yL = 0.05

    for _ in range(100):
        old_yH, old_yL = yH, yL

        # Low closed FWH balance:
        # (1-yH-yL)(h5-h4) = yL(h15-h7)
        delta_fw_low = stt[5]["h"] - stt[4]["h"]
        delta_steam_low = stt[15]["h"] - stt[7]["h"]
        yL = (1 - yH) * delta_fw_low / (delta_steam_low + delta_fw_low)

        # State 8 mixture:
        # mass entering from state 6 = 1-yH-yL
        # mass entering from state 7 = yL
        # mass leaving state 8 = 1-yH
        h8 = ((1 - yH - yL) * stt[6]["h"] + yL * stt[7]["h"]) / (1 - yH)

        stt[8] = {"P_kPa": Ph, "T_C": T_sat(Ph), "h": h8, "s": None, "v": None}

        # State 9 leaves high closed FWH
        stt[9] = {"P_kPa": Ph, "T_C": T_sat(Ph), "h": hf(Ph), "s": sf(Ph), "v": vf(Ph)}

        # High closed FWH balance:
        # (1-yH)(h9-h8) = yH(h14-h11)
        h11 = hf(Ph)
        delta_fw_high = stt[9]["h"] - stt[8]["h"]
        delta_steam_high = stt[14]["h"] - h11
        yH = delta_fw_high / (delta_steam_high + delta_fw_high)

        if abs(yH - old_yH) < 1e-8 and abs(yL - old_yL) < 1e-8:
            break

    # Recompute state 8 and 9 using final yH, yL
    h8 = ((1 - yH - yL) * stt[6]["h"] + yL * stt[7]["h"]) / (1 - yH)
    stt[8] = {"P_kPa": Ph, "T_C": T_sat(Ph), "h": h8, "s": None, "v": None}
    stt[9] = {"P_kPa": Ph, "T_C": T_sat(Ph), "h": hf(Ph), "s": sf(Ph), "v": vf(Ph)}

    # Pump main stream to boiler pressure
    stt[10] = {"P_kPa": Pb, "T_C": T_sat(Pb), "h": pump(stt[9]["h"], stt[9]["v"], Ph, Pb), "s": None, "v": stt[9]["v"]}

    # State 11: saturated liquid leaving high closed FWH
    stt[11] = {"P_kPa": Ph, "T_C": T_sat(Ph), "h": hf(Ph), "s": sf(Ph), "v": vf(Ph)}

    # Pump state 11 to boiler pressure before mixing
    h11_pumped = pump(stt[11]["h"], stt[11]["v"], Ph, Pb)

    # State 12 mixture:
    # mass from 10 = 1-yH
    # mass from pumped 11 = yH
    # total into boiler = 1
    h12 = (1 - yH) * stt[10]["h"] + yH * h11_pumped
    stt[12] = {"P_kPa": Pb, "T_C": T_sat(Pb), "h": h12, "s": None, "v": None}

    # Open FWH balance
    yO = ((1 - yH - yL) * (stt[3]["h"] - stt[2]["h"])) / (stt[18]["h"] - stt[2]["h"])

    # Turbine work pieces
    Wt_13_14 = stt[13]["h"] - stt[14]["h"]
    Wt_14_15 = (1 - yH) * (stt[14]["h"] - stt[15]["h"])
    Wt_15_16 = (1 - yH - yL) * (stt[15]["h"] - stt[16]["h"])
    Wt_17_18 = (1 - yH - yL) * (stt[17]["h"] - stt[18]["h"])
    Wt_18_19 = (1 - yH - yL - yO) * (stt[18]["h"] - stt[19]["h"])
    Wt_total = Wt_13_14 + Wt_14_15 + Wt_15_16 + Wt_17_18 + Wt_18_19

    # Pump work pieces
    Wp_1_2 = (1 - yH - yL - yO) * (stt[2]["h"] - stt[1]["h"])
    Wp_3_4 = (1 - yH - yL) * (stt[4]["h"] - stt[3]["h"])
    Wp_5_6 = (1 - yH - yL) * (stt[6]["h"] - stt[5]["h"])
    Wp_9_10 = (1 - yH) * (stt[10]["h"] - stt[9]["h"])
    Wp_11_12 = yH * (h11_pumped - stt[11]["h"])
    Wp_total = Wp_1_2 + Wp_3_4 + Wp_5_6 + Wp_9_10 + Wp_11_12

    # Heat transfer
    Q_boiler = stt[13]["h"] - stt[12]["h"]
    Q_reheat = (1 - yH - yL) * (stt[17]["h"] - stt[16]["h"])
    Q_in = Q_boiler + Q_reheat
    Q_cond = (1 - yH - yL - yO) * (stt[19]["h"] - stt[1]["h"])

    Wnet = Wt_total - Wp_total
    eta = Wnet / Q_in * 100

    # ---------------------------------------------------------
    # Mass flow rates by state
    # ---------------------------------------------------------
    # The cycle is normalized to 1 kg/s through the boiler.
    # Multiplying each fraction by m_dot_boiler gives actual kg/s.
    # For extraction states, the listed mass flow is the extracted stream.
    m_frac_by_state = {
        1: 1 - yH - yL - yO,   # condenser outlet
        2: 1 - yH - yL - yO,   # first pump outlet, entering open FWH
        3: 1 - yH - yL,        # open FWH outlet
        4: 1 - yH - yL,        # pump outlet to low closed FWH
        5: 1 - yH - yL,        # low closed FWH feedwater outlet
        6: 1 - yH - yL,        # pump outlet before low-drain mixing
        7: yL,                 # low closed FWH drain stream
        8: 1 - yH,             # mixture after low drain is added
        9: 1 - yH,             # high closed FWH feedwater outlet
        10: 1 - yH,            # pump outlet to boiler-pressure mixer
        11: yH,                # high closed FWH drain stream
        12: 1.0,               # mixture entering boiler
        13: 1.0,               # boiler outlet / turbine inlet
        14: yH,                # high-pressure extraction steam to high closed FWH
        15: yL,                # low-pressure extraction steam to low closed FWH
        16: 1 - yH - yL,       # flow to reheater
        17: 1 - yH - yL,       # reheater outlet
        18: yO,                # extraction steam to open FWH
        19: 1 - yH - yL - yO   # final turbine exhaust to condenser
    }

    state_table = pd.DataFrame.from_dict(stt, orient="index")
    state_table.index.name = "State"
    state_table = state_table.reset_index()
    state_table["m_dot_fraction"] = state_table["State"].map(m_frac_by_state)
    state_table["m_dot_kg_s"] = state_table["m_dot_fraction"] * m_dot_boiler
    state_table = state_table.sort_values("State")
    state_table = state_table.round(4)

    mass_flows = {
        "High closed FWH bleed fraction yH": yH,
        "Low closed FWH bleed fraction yL": yL,
        "Open FWH bleed fraction yO": yO,
        "Boiler mass flow rate (kg/s)": m_dot_boiler,
        "High closed FWH bleed flow, state 14 (kg/s)": yH * m_dot_boiler,
        "Low closed FWH bleed flow, state 15 (kg/s)": yL * m_dot_boiler,
        "Open FWH bleed flow, state 18 (kg/s)": yO * m_dot_boiler,
        "Flow after high bleed, states 8-10 (kg/s)": (1 - yH) * m_dot_boiler,
        "Flow after low bleed, states 3-6 and 16-17 (kg/s)": (1 - yH - yL) * m_dot_boiler,
        "Flow to condenser, states 1-2 and 19 (kg/s)": (1 - yH - yL - yO) * m_dot_boiler
    }

    component_results = {
        "Turbine Work 13→14 (kJ/kg)": Wt_13_14,
        "Turbine Work 14→15 (kJ/kg)": Wt_14_15,
        "Turbine Work 15→16 (kJ/kg)": Wt_15_16,
        "Turbine Work 17→18 (kJ/kg)": Wt_17_18,
        "Turbine Work 18→19 (kJ/kg)": Wt_18_19,
        "Total Turbine Work (kJ/kg)": Wt_total,
        "Pump Work 1→2 (kJ/kg)": Wp_1_2,
        "Pump Work 3→4 (kJ/kg)": Wp_3_4,
        "Pump Work 5→6 (kJ/kg)": Wp_5_6,
        "Pump Work 9→10 (kJ/kg)": Wp_9_10,
        "Pump Work 11→12 portion (kJ/kg)": Wp_11_12,
        "Total Pump Work (kJ/kg)": Wp_total,
        "Net Work (kJ/kg)": Wnet,
        "Boiler Heat Input (kJ/kg)": Q_boiler,
        "Reheater Heat Input (kJ/kg)": Q_reheat,
        "Total Heat Input (kJ/kg)": Q_in,
        "Condenser Heat Rejection (kJ/kg)": Q_cond,
        "Thermal Efficiency (%)": eta
    }

    return state_table, mass_flows, component_results


# -----------------------------
# Streamlit UI
# -----------------------------
st.set_page_config(page_title="Rankine Cycle Solver", layout="wide")
st.title("Rankine Reheat-Regenerative Cycle Solver")

with st.sidebar:
    st.header("Cycle Inputs")
    Pb = st.slider("Boiler Pressure (kPa)", 7000, 12000, 8000, 100)
    Pc = st.number_input("Condenser Pressure (kPa)", value=10)
    Pr = st.number_input("Reheat Pressure (kPa)", value=300)
    Ph = st.number_input("High Closed FWH Pressure (kPa)", value=6000)
    Pl = st.number_input("Low Closed FWH Pressure (kPa)", value=3500)
    Po = st.number_input("Open FWH Pressure (kPa)", value=100)
    Tin = st.number_input("Turbine Inlet Temperature (°C)", value=400)
    m_dot_boiler = st.number_input("Boiler Mass Flow Rate (kg/s)", value=1.0, min_value=0.0)

state_table, mass_flows, component_results = solve(Pb, Pc, Pr, Ph, Pl, Po, Tin, m_dot_boiler)

col1, col2, col3 = st.columns(3)
col1.metric("Thermal Efficiency", f"{component_results['Thermal Efficiency (%)']:.3f}%")
col2.metric("Net Work", f"{component_results['Net Work (kJ/kg)']:.2f} kJ/kg")
col3.metric("Heat Input", f"{component_results['Total Heat Input (kJ/kg)']:.2f} kJ/kg")

st.title("Rankine Reheat-Regenerative Cycle Solver")

st.subheader("Cycle Diagram")
st.image("cycle_diagram.png", caption="Reheat-Regenerative Rankine Cycle", width="stretch")

st.subheader("State Table")
st.dataframe(state_table, height=500, width="stretch")

st.subheader("Mass Flow Rates Summary")
mass_table = pd.DataFrame(mass_flows.items(), columns=["Quantity", "Value"]).round(6)
st.dataframe(mass_table, height=280, width="stretch")

st.subheader("Heat and Work Rate Calculations")
component_table = pd.DataFrame(component_results.items(), columns=["Quantity", "Value"]).round(5)
st.dataframe(component_table, height=500, width="stretch")

st.subheader("Efficiency vs Boiler Pressure")

Pb_range = np.linspace(7000, 12000, 25)
eff_boiler = []

for P in Pb_range:
    try:
        _, _, r = solve(Pb=P, Pc=Pc, Pr=Pr, Ph=Ph, Pl=Pl, Po=Po, Tin=Tin, m_dot_boiler=m_dot_boiler)
        eff_boiler.append(r["Thermal Efficiency (%)"])
    except:
        eff_boiler.append(np.nan)

fig1, ax1 = plt.subplots()
ax1.plot(Pb_range, eff_boiler, marker="o")
ax1.set_xlabel("Boiler Pressure (kPa)")
ax1.set_ylabel("Thermal Efficiency (%)")
ax1.grid(True)
st.pyplot(fig1)

boiler_table = pd.DataFrame({
    "Boiler Pressure (kPa)": Pb_range,
    "Efficiency (%)": eff_boiler
}).round(4)

st.dataframe(boiler_table, height=350, width="stretch")

st.subheader("Efficiency vs Open FWH Pressure")

Po_range = np.linspace(20, 300, 25)
eff_open = []

for Popen in Po_range:
    try:
        _, _, r = solve(Pb=Pb, Pc=Pc, Pr=Pr, Ph=Ph, Pl=Pl, Po=Popen, Tin=Tin, m_dot_boiler=m_dot_boiler)
        eff_open.append(r["Thermal Efficiency (%)"])
    except:
        eff_open.append(np.nan)
st.subheader("T-s Diagram")

# Use the current state table
ts_df = state_table.copy()

# Drop states with missing entropy or temperature
ts_plot = ts_df.dropna(subset=["s", "T_C"])

fig_ts, ax_ts = plt.subplots()

# Saturation dome
T_vals = np.linspace(5, 370, 300)
s_f = []
s_g = []

for T in T_vals:
    try:
        sf_val = CP.PropsSI("S", "T", T + 273.15, "Q", 0, WATER) / 1000
        sg_val = CP.PropsSI("S", "T", T + 273.15, "Q", 1, WATER) / 1000
        s_f.append(sf_val)
        s_g.append(sg_val)
    except:
        s_f.append(np.nan)
        s_g.append(np.nan)

ax_ts.plot(s_f, T_vals, label="Saturated Liquid Line")
ax_ts.plot(s_g, T_vals, label="Saturated Vapor Line")

# Cycle path
ax_ts.plot(
    ts_plot["s"],
    ts_plot["T_C"],
    marker="o",
    linewidth=2,
    label="Rankine Cycle"
)

# Label states
for _, row in ts_plot.iterrows():
    ax_ts.annotate(
        int(row["State"]),
        (row["s"], row["T_C"]),
        textcoords="offset points",
        xytext=(5, 5),
        fontsize=9
    )

ax_ts.set_xlabel("Entropy, s (kJ/kg·K)")
ax_ts.set_ylabel("Temperature, T (°C)")
ax_ts.set_title("T-s Diagram for Reheat-Regenerative Rankine Cycle")
ax_ts.grid(True)
ax_ts.legend()

st.pyplot(fig_ts)
fig2, ax2 = plt.subplots()
ax2.plot(Po_range, eff_open, marker="o")
ax2.set_xlabel("Open FWH Pressure (kPa)")
ax2.set_ylabel("Thermal Efficiency (%)")
ax2.grid(True)
st.pyplot(fig2)

open_table = pd.DataFrame({
    "Open FWH Pressure (kPa)": Po_range,
    "Efficiency (%)": eff_open
}).round(4)

st.dataframe(open_table, height=350, width="stretch")

best_index = np.nanargmax(eff_open)
st.success(
    f"Best open FWH pressure in this range: {Po_range[best_index]:.2f} kPa "
    f"with efficiency {eff_open[best_index]:.3f}%"
)