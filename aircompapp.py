# streamlit_app.py  (fixed with password protection using secrets)
import math
import pandas as pd
import streamlit as st
import CoolProp.CoolProp as CP

BAR_TO_PA = 1e5
PA_TO_BAR = 1e-5
C_TO_K = 273.15
SEC_PER_HOUR = 3600.0

st.set_page_config(page_title="Multi-Stage Air Compressors (O₂/N₂ + CoolProp)", layout="wide")

# -----------------------------
# Password Protection
# -----------------------------
def check_password():
    """Returns `True` if the user has entered the correct password."""
    
    def password_entered():
        """Checks whether a password entered by the user is correct."""
        if st.session_state["password"] == st.secrets["password"]:
            st.session_state["password_correct"] = True
            del st.session_state["password"]  # Don't store password
        else:
            st.session_state["password_correct"] = False

    # First run or password not correct
    if "password_correct" not in st.session_state:
        # Show input for password
        st.text_input(
            "Enter Password", type="password", on_change=password_entered, key="password"
        )
        st.info("Please enter the password to access the application.")
        return False
    elif not st.session_state["password_correct"]:
        # Password incorrect
        st.text_input(
            "Enter Password", type="password", on_change=password_entered, key="password"
        )
        st.error("😕 Password incorrect. Please try again.")
        return False
    else:
        # Password correct
        return True

if not check_password():
    st.stop()

# -----------------------------
# Main Application (after password check)
# -----------------------------

st.title("Multi-Stage Air Compressor Calculator — O₂/N₂ (CoolProp)")
st.caption("Real-gas properties via CoolProp (HEOS). Intercooling to 30 °C. Power via isentropic compression with user-selected efficiency.")

# -----------------------------
# Helpers
# -----------------------------
def make_mixture(o2_mole_frac: float):
    """Create an O2/N2 CoolProp mixture state object using HEOS backend."""
    xO2 = max(0.0, min(1.0, o2_mole_frac))
    xN2 = 1.0 - xO2
    asys = CP.AbstractState("HEOS", "Nitrogen&Oxygen")
    # Order must match the components string
    asys.set_mole_fractions([xN2, xO2])
    return asys

def state_PT(asys, P_bar: float, T_C: float):
    P = P_bar * BAR_TO_PA
    T = T_C + C_TO_K
    asys.update(CP.PT_INPUTS, P, T)
    return {
        "P_bar": P_bar,
        "T_C": T_C,
        "rho": asys.rhomass(),          # kg/m^3
        "h": asys.hmass(),              # J/kg
        "s": asys.smass(),              # J/kg/K
        "T_K": T,
        "P_Pa": P,
    }

def rho_PT(asys, P_bar: float, T_C: float) -> float:
    P = P_bar * BAR_TO_PA
    T = T_C + C_TO_K
    asys.update(CP.PT_INPUTS, P, T)
    return asys.rhomass()

def isentropic_h2(asys, P2_bar: float, s_target: float) -> float:
    """
    Robustly find the isentropic target enthalpy h2s at pressure P2 for a mixture
    by solving s(P2, T) = s_target using PT inputs (avoids PSmass flash).
    """
    P2 = P2_bar * BAR_TO_PA

    def s_at_T(TK: float) -> float:
        asys.update(CP.PT_INPUTS, P2, TK)
        return asys.smass()

    # Start near the current temperature if the state was just updated
    T0 = asys.T() if asys.T() > 1.0 else (300.0)  # K

    # Build a bracket [T_low, T_high] where s crosses s_target
    T_low = max(180.0, 0.5 * T0)      # avoid cryogenic two-phase region
    T_high = min(3000.0, 2.5 * T0)    # keep within sane bounds

    # Expand the bracket if needed
    max_expand = 60
    expanded = 0
    while True:
        try:
            s_low = s_at_T(T_low)
            s_high = s_at_T(T_high)
            if (s_low - s_target) * (s_high - s_target) <= 0.0:
                break  # bracket found
            # Expand outward
            T_low = max(120.0, T_low * 0.8)
            T_high = min(4000.0, T_high * 1.25)
            expanded += 1
            if expanded > max_expand:
                # Last resort: nudge around T0
                T_low = max(200.0, T0 * 0.7)
                T_high = min(3500.0, T0 * 1.6)
                s_low = s_at_T(T_low)
                s_high = s_at_T(T_high)
                break
        except Exception:
            # If we hit a bad state, gently move away
            T_low = max(150.0, T_low + 10.0)
            T_high = min(3500.0, T_high + 20.0)
            expanded += 1
            if expanded > max_expand:
                break

    # Bisection to solve s(P2,T) = s_target
    iter_max = 80
    for _ in range(iter_max):
        T_mid = 0.5 * (T_low + T_high)
        s_mid = s_at_T(T_mid)
        if abs(s_mid - s_target) < 1e-6 or (T_high - T_low) < 1e-6:
            T_star = T_mid
            break
        # Keep the sub-interval that brackets the root
        if (s_low - s_target) * (s_mid - s_target) <= 0.0:
            T_high = T_mid
            s_high = s_mid
        else:
            T_low = T_mid
            s_low = s_mid
    else:
        T_star = T_mid

    # Return h(P2, T_star)
    asys.update(CP.PT_INPUTS, P2, T_star)
    return asys.hmass()

def T_from_HP(asys, h: float, P_bar: float) -> float:
    P = P_bar * BAR_TO_PA
    asys.update(CP.HmassP_INPUTS, h, P)
    return asys.T() - C_TO_K  # °C

def compress_multistage(
    asys,
    m_dot: float,
    P_in_bar: float,
    T_in_C: float,
    P_out_bar: float,
    n_stages: int,
    eta_s: float,
    intercool_T_C: float,
):
    """
    Multi-stage compressor with equal pressure ratio per stage and intercooling to intercool_T_C after each stage
    (including the final discharge cooler).
    Returns totals and a per-stage DataFrame.
    """
    n_stages = max(1, int(n_stages))
    eta_s = max(1e-6, min(1.0, eta_s))
    r_total = P_out_bar / P_in_bar
    r_stage = r_total ** (1.0 / n_stages)

    stage_rows = []
    total_power_W = 0.0

    P1_bar = P_in_bar
    T1_C = T_in_C
    # Precompute inlet state
    st1 = state_PT(asys, P1_bar, T1_C)
    h1 = st1["h"]
    s1 = st1["s"]

    for i in range(1, n_stages + 1):
        P2_bar = P1_bar * r_stage

        # Isentropic target enthalpy at P2 via PT rootfind
        h2s = isentropic_h2(asys, P2_bar, s1)
        # Actual enthalpy rise with isentropic efficiency
        h2 = h1 + (h2s - h1) / eta_s
        # Actual discharge temperature
        T2_C = T_from_HP(asys, h2, P2_bar)

        # Stage power
        Wdot_stage_W = m_dot * (h2 - h1)  # W
        total_power_W += Wdot_stage_W

        # Volumetric flows (discharge actual, then aftercool to intercool_T_C)
        rho2 = rho_PT(asys, P2_bar, T2_C)
        Vdot2_m3s = m_dot / max(1e-12, rho2)
        Vdot2_m3h = Vdot2_m3s * SEC_PER_HOUR

        # Aftercool to intercool_T_C at same pressure
        rho_ac = rho_PT(asys, P2_bar, intercool_T_C)
        Vdot_ac_m3s = m_dot / max(1e-12, rho_ac)
        Vdot_ac_m3h = Vdot_ac_m3s * SEC_PER_HOUR

        stage_rows.append(
            {
                "Stage": i,
                "P_in (bar)": P1_bar,
                "T_in (°C)": T1_C,
                "P_out (bar)": P2_bar,
                "T_out (°C)": T2_C,
                "Power (kW)": Wdot_stage_W / 1e3,
                "V̇_out @T_out (m³/h)": Vdot2_m3h,
                "V̇_aftercool @30 °C (m³/h)": Vdot_ac_m3h,
            }
        )

        # Intercool to next inlet
        P1_bar = P2_bar
        T1_C = intercool_T_C
        st1 = state_PT(asys, P1_bar, T1_C)
        h1 = st1["h"]
        s1 = st1["s"]

    # Final discharge (before final cooler) is last stage's T_out
    final_discharge_T_C = stage_rows[-1]["T_out (°C)"]
    # Final aftercool volumetric flow (at P_out_bar, intercool_T_C)
    rho_final_cool = rho_PT(asys, P_out_bar, intercool_T_C)
    Vdot_final_cool_m3h = (m_dot / max(1e-12, rho_final_cool)) * SEC_PER_HOUR

    df = pd.DataFrame(stage_rows)
    return {
        "total_power_MW": total_power_W / 1e6,
        "final_discharge_T_C": final_discharge_T_C,
        "Vdot_aftercool_m3h": Vdot_final_cool_m3h,
        "stages_df": df,
    }

# -----------------------------
# Sidebar controls (mixture & method)
# -----------------------------
st.sidebar.header("Mixture & Method")
o2_frac = st.sidebar.number_input("O₂ mole fraction", min_value=0.0, max_value=1.0, value=0.2095, step=0.0005, format="%.4f")
eta_s = st.sidebar.slider("Isentropic efficiency (compressor)", min_value=0.5, max_value=1.0, value=0.80, step=0.01)
T_intercool_C = st.sidebar.number_input("Intercooler outlet temperature (°C)", value=30.0, step=1.0)

asys = make_mixture(o2_frac)

# -----------------------------
# Compressor 1: 1 bar, 20 °C -> 7 bar (default), 2 stages
# -----------------------------
st.subheader("Compressor 1 — Inlet: 20 °C & 1 bar → Outlet: default 7 bar (changeable)")
col1, col2, col3, col4 = st.columns(4)
with col1:
    Vdot1_in_m3h = st.number_input("Inlet volume flow (m³/h) @ 20 °C, 1 bar", min_value=0.0, value=10000000.0, step=10.0, format="%.3f")
with col2:
    P1_in_bar = st.number_input("Inlet pressure (bar)", min_value=0.1, value=1.0, step=0.1, format="%.3f")
with col3:
    T1_in_C = st.number_input("Inlet temperature (°C)", value=20.0, step=0.5, format="%.2f")
with col4:
    P1_out_bar = st.number_input("Outlet pressure (bar)", min_value=0.2, value=7.0, step=0.1, format="%.3f")

stages_1 = st.number_input("Number of stages (Compressor 1)", min_value=1, value=2, step=1)

# Mass flow from inlet volumetric flow
rho_in_1 = rho_PT(asys, P1_in_bar, T1_in_C)
Vdot1_in_m3s = Vdot1_in_m3h / SEC_PER_HOUR
m_dot_1 = rho_in_1 * Vdot1_in_m3s  # kg/s

res1 = compress_multistage(
    asys=asys,
    m_dot=m_dot_1,
    P_in_bar=P1_in_bar,
    T_in_C=T1_in_C,
    P_out_bar=P1_out_bar,
    n_stages=stages_1,
    eta_s=eta_s,
    intercool_T_C=T_intercool_C,
)

# Display summary for Compressor 1
mcol1, mcol2, mcol3 = st.columns(3)
mcol1.metric("Total shaft power (MW)", f"{res1['total_power_MW']:.6f}")
mcol2.metric("Final discharge temperature (°C)", f"{res1['final_discharge_T_C']:.2f}")
mcol3.metric(f"Outlet volume aftercool @ {T_intercool_C:.0f} °C (m³/h)", f"{res1['Vdot_aftercool_m3h']:.3f}")

st.dataframe(res1["stages_df"].style.format({
    "P_in (bar)": "{:.4f}",
    "T_in (°C)": "{:.2f}",
    "P_out (bar)": "{:.4f}",
    "T_out (°C)": "{:.2f}",
    "Power (kW)": "{:.3f}",
    "V̇_out @T_out (m³/h)": "{:.3f}",
    "V̇_aftercool @30 °C (m³/h)": "{:.3f}",
}), use_container_width=True)

st.divider()

# -----------------------------
# Compressor 2: feed from 7 bar, 30 °C (auto) -> user pressure (default 45 bar), default 3 stages
# -----------------------------
st.subheader("Compressor 2 — Inlet: use 7 bar stream cooled to 30 °C → Outlet: default 45 bar (changeable)")

default_Vdot2_in = float(res1["Vdot_aftercool_m3h"]) if math.isfinite(res1["Vdot_aftercool_m3h"]) else 0.0

c2col1, c2col2, c2col3, c2col4 = st.columns(4)
with c2col1:
    Vdot2_in_m3h = st.number_input("Inlet volume flow (m³/h) @ 30 °C, 7 bar", min_value=0.0, value=round(default_Vdot2_in, 6), step=1.0, format="%.6f")
with c2col2:
    P2_in_bar = st.number_input("Inlet pressure (bar)", min_value=0.1, value=float(P1_out_bar), step=0.1, format="%.3f")
with c2col3:
    T2_in_C = st.number_input("Inlet temperature (°C)", value=float(T_intercool_C), step=0.5, format="%.2f")
with c2col4:
    P2_out_bar = st.number_input("Outlet pressure (bar)", min_value=0.2, value=45.0, step=0.5, format="%.3f")

stages_2 = st.number_input("Number of stages (Compressor 2)", min_value=1, value=3, step=1)

# Mass flow for compressor 2 based on its inlet state and inlet volumetric flow
rho_in_2 = rho_PT(asys, P2_in_bar, T2_in_C)
Vdot2_in_m3s = Vdot2_in_m3h / SEC_PER_HOUR
m_dot_2 = rho_in_2 * Vdot2_in_m3s  # kg/s

res2 = compress_multistage(
    asys=asys,
    m_dot=m_dot_2,
    P_in_bar=P2_in_bar,
    T_in_C=T2_in_C,
    P_out_bar=P2_out_bar,
    n_stages=stages_2,
    eta_s=eta_s,
    intercool_T_C=T_intercool_C,
)

# Display summary for Compressor 2
ncol1, ncol2, ncol3 = st.columns(3)
ncol1.metric("Total shaft power (MW)", f"{res2['total_power_MW']:.6f}")
ncol2.metric("Final discharge temperature (°C)", f"{res2['final_discharge_T_C']:.2f}")
ncol3.metric(f"Outlet volume aftercool @ {T_intercool_C:.0f} °C (m³/h)", f"{res2['Vdot_aftercool_m3h']:.6f}")

st.dataframe(res2["stages_df"].style.format({
    "P_in (bar)": "{:.4f}",
    "T_in (°C)": "{:.2f}",
    "P_out (bar)": "{:.4f}",
    "T_out (°C)": "{:.2f}",
    "Power (kW)": "{:.6f}",
    "V̇_out @T_out (m³/h)": "{:.6f}",
    "V̇_aftercool @30 °C (m³/h)": "{:.6f}",
}), use_container_width=True)

# -----------------------------
# Notes / assumptions
# -----------------------------
with st.expander("Modeling notes & assumptions"):
    st.markdown("""
- Composition is **O₂ + N₂ only** (CoolProp HEOS binary mixture). Default O₂ mole fraction is **0.2095** (dry air).
- Each compressor uses **equal pressure ratio per stage** with **isentropic compression** and a user-selected **isentropic efficiency** (default 0.80).
- **Intercooler / aftercooler** returns the gas to **set temperature** (default **30 °C**) at the respective stage outlet pressure.  
- Power is **shaft power** from enthalpy rise:  \n
  \\( \\dot W = \\dot m\\,(h_{2}-h_{1}) \\) with \\( h_{2} = h_{1} + (h_{2s}-h_{1})/\\eta_s \\).  
- Isentropic target is found by solving **s(P₂,T₂)=s₁** with a numerical bracket + bisection (no `PSmass` flash needed).
- Volumetric flows use \\(\\dot V = \\dot m/\\rho(P,T)\\) with real-gas density from CoolProp.
- Units: pressures in **bar**, temperatures in **°C**, power in **MW**, volumes in **m³/h**.
""")