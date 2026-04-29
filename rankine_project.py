import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

WATER = "Water"

def h_PT(P, T):
    return CP.PropsSI("H", "P", P*1000, "T", T+273.15, WATER)/1000

def s_PT(P, T):
    return CP.PropsSI("S", "P", P*1000, "T", T+273.15, WATER)/1000

def h_Ps(P, s):
    return CP.PropsSI("H", "P", P*1000, "S", s*1000, WATER)/1000

def hf(P):
    return CP.PropsSI("H", "P", P*1000, "Q", 0, WATER)/1000

def sf(P):
    return CP.PropsSI("S", "P", P*1000, "Q", 0, WATER)/1000

def vf(P):
    rho = CP.PropsSI("D", "P", P*1000, "Q", 0, WATER)
    return 1/rho

def T_sat(P):
    return CP.PropsSI("T", "P", P*1000, "Q", 0, WATER) - 273.15

def T_Ps(P, s):
    return CP.PropsSI("T", "P", P*1000, "S", s*1000, WATER) - 273.15

def pump(h, v, P1, P2):
    return h + v*(P2-P1)

def solve(Pb=8000, Pc=10, Pr=300, Ph=6000, Pl=3500, Po=100, Tin=400):
    st = {}

    st[1] = {"P_kPa": Pc, "T_C": T_sat(Pc), "h": hf(Pc), "s": sf(Pc), "v": vf(Pc)}
    st[2] = {"P_kPa": Po, "T_C": T_sat(Po), "h": pump(st[1]["h"], st[1]["v"], Pc, Po), "s": None, "v": st[1]["v"]}

    st[3] = {"P_kPa": Po, "T_C": T_sat(Po), "h": hf(Po), "s": sf(Po), "v": vf(Po)}
    st[4] = {"P_kPa": Pl, "T_C": T_sat(Pl), "h": pump(st[3]["h"], st[3]["v"], Po, Pl), "s": None, "v": st[3]["v"]}

    st[5] = {"P_kPa": Pl, "T_C": T_sat(Pl), "h": hf(Pl), "s": sf(Pl), "v": vf(Pl)}
    st[6] = {"P_kPa": Ph, "T_C": T_sat(Ph), "h": pump(st[5]["h"], st[5]["v"], Pl, Ph), "s": None, "v": st[5]["v"]}

    st[7] = {"P_kPa": Pl, "T_C": T_sat(Pl), "h": hf(Pl), "s": sf(Pl), "v": vf(Pl)}

    st[13] = {"P_kPa": Pb, "T_C": Tin, "h": h_PT(Pb, Tin), "s": s_PT(Pb, Tin), "v": None}
    st[14] = {"P_kPa": Ph, "T_C": T_Ps(Ph, st[13]["s"]), "h": h_Ps(Ph, st[13]["s"]), "s": st[13]["s"], "v": None}
    st[15] = {"P_kPa": Pl, "T_C": T_Ps(Pl, st[13]["s"]), "h": h_Ps(Pl, st[13]["s"]), "s": st[13]["s"], "v": None}
    st[16] = {"P_kPa": Pr, "T_C": T_Ps(Pr, st[13]["s"]), "h": h_Ps(Pr, st[13]["s"]), "s": st[13]["s"], "v": None}

    st[17] = {"P_kPa": Pr, "T_C": Tin, "h": h_PT(Pr, Tin), "s": s_PT(Pr, Tin), "v": None}
    st[18] = {"P_kPa": Po, "T_C": T_Ps(Po, st[17]["s"]), "h": h_Ps(Po, st[17]["s"]), "s": st[17]["s"], "v": None}
    st[19] = {"P_kPa": Pc, "T_C": T_Ps(Pc, st[17]["s"]), "h": h_Ps(Pc, st[17]["s"]), "s": st[17]["s"], "v": None}

    yL = (st[5]["h"] - st[4]["h"]) / (st[15]["h"] - st[7]["h"])

    st[8] = {
        "P_kPa": Ph,
        "T_C": T_sat(Ph),
        "h": (st[6]["h"] + yL*st[7]["h"]) / (1+yL),
        "s": None,
        "v": None
    }

    st[9] = {"P_kPa": Ph, "T_C": T_sat(Ph), "h": hf(Ph), "s": sf(Ph), "v": vf(Ph)}

    yH = (st[9]["h"] - st[8]["h"]) / (st[14]["h"] - hf(Ph))

    st[10] = {
        "P_kPa": Pb,
        "T_C": T_sat(Pb),
        "h": pump(st[9]["h"], st[9]["v"], Ph, Pb),
        "s": None,
        "v": st[9]["v"]
    }

    st[11] = {"P_kPa": Ph, "T_C": T_sat(Ph), "h": hf(Ph), "s": sf(Ph), "v": vf(Ph)}

    h11_pumped = pump(st[11]["h"], st[11]["v"], Ph, Pb)

    st[12] = {
        "P_kPa": Pb,
        "T_C": T_sat(Pb),
        "h": (st[10]["h"] + yH*h11_pumped) / (1+yH),
        "s": None,
        "v": None
    }

    yO = ((1-yH-yL)*(st[3]["h"] - st[2]["h"])) / (st[18]["h"] - st[2]["h"])

    Wt_13_14 = st[13]["h"] - st[14]["h"]
    Wt_14_15 = (1-yH)*(st[14]["h"] - st[15]["h"])
    Wt_15_16 = (1-yH-yL)*(st[15]["h"] - st[16]["h"])
    Wt_17_18 = (1-yH-yL)*(st[17]["h"] - st[18]["h"])
    Wt_18_19 = (1-yH-yL-yO)*(st[18]["h"] - st[19]["h"])

    Wt = Wt_13_14 + Wt_14_15 + Wt_15_16 + Wt_17_18 + Wt_18_19

    Wp_1_2 = (1-yH-yL-yO)*(st[2]["h"] - st[1]["h"])
    Wp_3_4 = (1-yH-yL)*(st[4]["h"] - st[3]["h"])
    Wp_5_6 = st[6]["h"] - st[5]["h"]
    Wp_9_10 = st[10]["h"] - st[9]["h"]
    Wp_11_12 = yH*(h11_pumped - st[11]["h"])

    Wp = Wp_1_2 + Wp_3_4 + Wp_5_6 + Wp_9_10 + Wp_11_12

    Q_boiler = st[13]["h"] - st[12]["h"]
    Q_reheat = (1-yH-yL)*(st[17]["h"] - st[16]["h"])
    Qin = Q_boiler + Q_reheat

    Q_condenser = (1-yH-yL-yO)*(st[19]["h"] - st[1]["h"])

    Wnet = Wt - Wp
    eta = Wnet / Qin * 100

    mass_flows = {
        "Boiler Mass Flow": 1.0,
        "High Closed FWH Bleed yH": yH,
        "Low Closed FWH Bleed yL": yL,
        "Open FWH Bleed yO": yO,
        "Flow After High Bleed": 1-yH,
        "Flow After Low Bleed": 1-yH-yL,
        "Flow To Condenser": 1-yH-yL-yO
    }

    component_results = {
        "Turbine Work 13→14 (kJ/kg)": Wt_13_14,
        "Turbine Work 14→15 (kJ/kg)": Wt_14_15,
        "Turbine Work 15→16 (kJ/kg)": Wt_15_16,
        "Turbine Work 17→18 (kJ/kg)": Wt_17_18,
        "Turbine Work 18→19 (kJ/kg)": Wt_18_19,
        "Total Turbine Work (kJ/kg)": Wt,

        "Pump Work 1→2 (kJ/kg)": Wp_1_2,
        "Pump Work 3→4 (kJ/kg)": Wp_3_4,
        "Pump Work 5→6 (kJ/kg)": Wp_5_6,
        "Pump Work 9→10 (kJ/kg)": Wp_9_10,
        "Pump Work 11→12 portion (kJ/kg)": Wp_11_12,
        "Total Pump Work (kJ/kg)": Wp,

        "Net Work (kJ/kg)": Wnet,
        "Boiler Heat Input (kJ/kg)": Q_boiler,
        "Reheater Heat Input (kJ/kg)": Q_reheat,
        "Total Heat Input (kJ/kg)": Qin,
        "Condenser Heat Rejection (kJ/kg)": Q_condenser,
        "Thermal Efficiency (%)": eta
    }

    state_table = pd.DataFrame.from_dict(st, orient="index")
    state_table.index.name = "State"
    state_table = state_table.reset_index()
    state_table = state_table.sort_values("State")
    state_table = state_table.round(4)

    return state_table, mass_flows, component_results


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

state_table, mass_flows, component_results = solve(Pb, Pc, Pr, Ph, Pl, Po, Tin)

col1, col2, col3 = st.columns(3)
col1.metric("Thermal Efficiency", f"{component_results['Thermal Efficiency (%)']:.3f}%")
col2.metric("Net Work", f"{component_results['Net Work (kJ/kg)']:.2f} kJ/kg")
col3.metric("Heat Input", f"{component_results['Total Heat Input (kJ/kg)']:.2f} kJ/kg")

st.subheader("State Table")
st.dataframe(state_table, height=500, width="stretch")

st.subheader("Mass Flowrates")
mass_table = pd.DataFrame(mass_flows.items(), columns=["Quantity", "Value"]).round(6)
st.dataframe(mass_table, height=280, width="stretch")

st.subheader("Heat and Work Rate Calculations")
component_table = pd.DataFrame(component_results.items(), columns=["Quantity", "Value"]).round(5)
st.dataframe(component_table, height=500, width="stretch")

st.subheader("Efficiency vs Boiler Pressure")

Pb_range = np.linspace(7000, 12000, 25)
eff = []

for P in Pb_range:
    try:
        _, _, r = solve(Pb=P, Pc=Pc, Pr=Pr, Ph=Ph, Pl=Pl, Po=Po, Tin=Tin)
        eff.append(r["Thermal Efficiency (%)"])
    except:
        eff.append(np.nan)

fig, ax = plt.subplots()
ax.plot(Pb_range, eff, marker="o")
ax.set_xlabel("Boiler Pressure (kPa)")
ax.set_ylabel("Thermal Efficiency (%)")
ax.grid(True)

st.pyplot(fig)
st.subheader("Boiler Pressure Study Table")
study_table = pd.DataFrame({
    "Boiler Pressure (kPa)": Pb_range,
    "Efficiency (%)": eff
}).round(4)

st.dataframe(study_table, height=350, width="stretch")