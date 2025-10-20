# streamlit_app.py - Multi-Stage Air Compressor with O₂/N₂ + CoolProp
import math
import pandas as pd
import streamlit as st
import CoolProp.CoolProp as CP

# Constants
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

st.title("Multi-Stage Air Compressor with Intercooling")

# Sidebar inputs
st.sidebar.header("Input Parameters")

# Gas composition (air: 21% O2, 79% N2 by volume)
o2_frac = st.sidebar.slider("O₂ Mole Fraction", 0.0, 1.0, 0.21, 0.01)
n2_frac = 1 - o2_frac
st.sidebar.write(f"N₂ Mole Fraction: {n2_frac:.2f}")

# Inlet conditions
vol_flow_inlet = st.sidebar.number_input("Inlet Volumetric Flow (m³/h)", 
                                          min_value=1.0, value=1000.0, step=10.0)
T_inlet = st.sidebar.number_input("Inlet Temperature (°C)", value=20.0)
P_inlet = st.sidebar.number_input("Inlet Pressure (bar)", value=1.0, step=0.1)

# First stage compression
st.sidebar.header("First Compression Section")
stages_1 = st.sidebar.number_input("Number of Stages (Section 1)", 
                                    min_value=1, value=2, step=1)
P_intermediate = st.sidebar.number_input("Intermediate Pressure (bar)", 
                                          min_value=1.1, value=7.0, step=0.5)
T_aftercool_1 = st.sidebar.number_input("After-cooler Temperature (°C)", value=30.0)

# Second stage compression
st.sidebar.header("Second Compression Section")
stages_2 = st.sidebar.number_input("Number of Stages (Section 2)", 
                                    min_value=1, value=3, step=1)
P_outlet = st.sidebar.number_input("Final Outlet Pressure (bar)", 
                                    min_value=1.1, value=45.0, step=1.0)

# Efficiency
st.sidebar.header("Efficiency Parameters")
eta_isentropic = st.sidebar.slider("Isentropic Efficiency", 0.5, 1.0, 0.85, 0.01)
eta_mech = st.sidebar.slider("Mechanical Efficiency", 0.9, 1.0, 0.98, 0.01)

st.sidebar.info(f"""
**Efficiency Explanation:**
- **Isentropic Efficiency** (η_isen = {eta_isentropic:.0%}): Accounts for thermodynamic losses in compression
- **Mechanical Efficiency** (η_mech = {eta_mech:.0%}): Accounts for mechanical losses (bearings, seals, etc.)
- **Shaft Power** = Gas Power / η_mech
""")

# Gas properties functions
def get_gas_properties(T_C, P_bar, x_o2, x_n2):
    """Calculate mixture properties using CoolProp"""
    T_K = T_C + C_TO_K
    P_Pa = P_bar * BAR_TO_PA
    
    # Molecular weights
    MW_O2 = 31.999  # g/mol
    MW_N2 = 28.014  # g/mol
    MW_mix = x_o2 * MW_O2 + x_n2 * MW_N2
    
    # Get properties for each component
    cp_O2 = CP.PropsSI('C', 'T', T_K, 'P', P_Pa, 'O2')
    cp_N2 = CP.PropsSI('C', 'T', T_K, 'P', P_Pa, 'N2')
    cp_mix = x_o2 * cp_O2 + x_n2 * cp_N2
    
    cv_O2 = CP.PropsSI('O', 'T', T_K, 'P', P_Pa, 'O2')
    cv_N2 = CP.PropsSI('O', 'T', T_K, 'P', P_Pa, 'N2')
    cv_mix = x_o2 * cv_O2 + x_n2 * cv_N2
    
    gamma = cp_mix / cv_mix
    
    # Density
    rho_O2 = CP.PropsSI('D', 'T', T_K, 'P', P_Pa, 'O2')
    rho_N2 = CP.PropsSI('D', 'T', T_K, 'P', P_Pa, 'N2')
    rho_mix = 1 / (x_o2/rho_O2 + x_n2/rho_N2)
    
    return MW_mix, cp_mix, gamma, rho_mix

def compress_stage(T_in_C, P_in_bar, P_out_bar, mass_flow, x_o2, x_n2, eta_isen, eta_mech):
    """Compress gas through one stage with both efficiencies"""
    T_in_K = T_in_C + C_TO_K
    
    MW, cp, gamma, rho_in = get_gas_properties(T_in_C, P_in_bar, x_o2, x_n2)
    
    # Isentropic temperature rise
    PR = P_out_bar / P_in_bar
    T_out_s_K = T_in_K * (PR ** ((gamma - 1) / gamma))
    
    # Actual temperature with isentropic efficiency
    T_out_K = T_in_K + (T_out_s_K - T_in_K) / eta_isen
    T_out_C = T_out_K - C_TO_K
    
    # Gas power (thermodynamic power transferred to gas)
    power_gas_W = mass_flow * cp * (T_out_K - T_in_K)
    
    # Shaft power (actual power required at shaft, accounting for mechanical losses)
    power_shaft_W = power_gas_W / eta_mech
    
    # Convert to MW
    power_gas_MW = power_gas_W / 1e6
    power_shaft_MW = power_shaft_W / 1e6
    
    return T_out_C, power_gas_MW, power_shaft_MW

def compress_multi_stage(T_in_C, P_in_bar, P_out_bar, vol_flow_m3h, 
                         n_stages, x_o2, x_n2, eta_isen, eta_mech):
    """Compress through multiple stages with equal pressure ratios"""
    # Calculate mass flow
    MW, cp, gamma, rho_in = get_gas_properties(T_in_C, P_in_bar, x_o2, x_n2)
    vol_flow_m3s = vol_flow_m3h / SEC_PER_HOUR
    mass_flow = rho_in * vol_flow_m3s  # kg/s
    
    # Pressure ratio per stage
    PR_total = P_out_bar / P_in_bar
    PR_stage = PR_total ** (1 / n_stages)
    
    T_current = T_in_C
    P_current = P_in_bar
    total_power_gas = 0
    total_power_shaft = 0
    
    stage_results = []
    
    for i in range(n_stages):
        P_next = P_current * PR_stage
        T_out, power_gas, power_shaft = compress_stage(
            T_current, P_current, P_next, 
            mass_flow, x_o2, x_n2, eta_isen, eta_mech
        )
        
        total_power_gas += power_gas
        total_power_shaft += power_shaft
        
        stage_results.append({
            'stage': i + 1,
            'P_in': P_current,
            'P_out': P_next,
            'T_in': T_current,
            'T_out': T_out,
            'power_gas': power_gas,
            'power_shaft': power_shaft
        })
        
        # Intercooling (back to inlet temp for intermediate stages)
        if i < n_stages - 1:
            T_current = T_in_C  # Intercooler
        else:
            T_current = T_out
        P_current = P_next
    
    return T_current, total_power_gas, total_power_shaft, mass_flow, stage_results

# Create tabs
tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Calculator", 
    "Table P=45 bar", 
    "Table P=50 bar", 
    "Table P=70 bar", 
    "Table P=84 bar"
])

with tab1:
    # Main calculations
    st.header("Compression Results")
    
    # Section 1: Compress to intermediate pressure
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("Section 1: Compression to Intermediate Pressure")
        
        T_out_1, power_gas_1, power_shaft_1, mass_flow, stages_data_1 = compress_multi_stage(
            T_inlet, P_inlet, P_intermediate, vol_flow_inlet, 
            stages_1, o2_frac, n2_frac, eta_isentropic, eta_mech
        )
        
        st.metric("Outlet Temperature", f"{T_out_1:.1f} °C")
        st.metric("Outlet Pressure", f"{P_intermediate:.1f} bar")
        st.metric("Gas Power (Thermodynamic)", f"{power_gas_1:.3f} MW")
        st.metric("Shaft Power Required", f"{power_shaft_1:.3f} MW", 
                  help="Includes mechanical losses")
        st.metric("Mass Flow Rate", f"{mass_flow:.2f} kg/s")
        
        # Calculate outlet volume flow before cooling
        _, _, _, rho_out_1 = get_gas_properties(T_out_1, P_intermediate, o2_frac, n2_frac)
        vol_flow_out_1 = (mass_flow / rho_out_1) * SEC_PER_HOUR
        st.metric("Outlet Vol Flow (before cooling)", f"{vol_flow_out_1:.1f} m³/h")
        
        # After cooling
        _, _, _, rho_cool_1 = get_gas_properties(T_aftercool_1, P_intermediate, o2_frac, n2_frac)
        vol_flow_cool_1 = (mass_flow / rho_cool_1) * SEC_PER_HOUR
        st.metric("Vol Flow after cooling to 30°C", f"{vol_flow_cool_1:.1f} m³/h")
        
        with st.expander("Stage Details - Section 1"):
            for stage in stages_data_1:
                st.write(f"**Stage {stage['stage']}**")
                st.write(f"  Pressure: {stage['P_in']:.2f} → {stage['P_out']:.2f} bar")
                st.write(f"  Temperature: {stage['T_in']:.1f} → {stage['T_out']:.1f} °C")
                st.write(f"  Gas Power: {stage['power_gas']:.3f} MW")
                st.write(f"  Shaft Power: {stage['power_shaft']:.3f} MW")
                st.write("---")
    
    with col2:
        st.subheader("Section 2: Compression to Final Pressure")
        
        # Section 2 uses cooled gas from section 1
        T_out_2, power_gas_2, power_shaft_2, mass_flow_2, stages_data_2 = compress_multi_stage(
            T_aftercool_1, P_intermediate, P_outlet, vol_flow_cool_1, 
            stages_2, o2_frac, n2_frac, eta_isentropic, eta_mech
        )
        
        st.metric("Outlet Temperature", f"{T_out_2:.1f} °C")
        st.metric("Outlet Pressure", f"{P_outlet:.1f} bar")
        st.metric("Gas Power (Thermodynamic)", f"{power_gas_2:.3f} MW")
        st.metric("Shaft Power Required", f"{power_shaft_2:.3f} MW",
                  help="Includes mechanical losses")
        st.metric("Mass Flow Rate", f"{mass_flow_2:.2f} kg/s")
        
        # Calculate outlet volume flow before final cooling
        _, _, _, rho_out_2 = get_gas_properties(T_out_2, P_outlet, o2_frac, n2_frac)
        vol_flow_out_2 = (mass_flow_2 / rho_out_2) * SEC_PER_HOUR
        st.metric("Outlet Vol Flow (before cooling)", f"{vol_flow_out_2:.1f} m³/h")
        
        # After final cooling
        _, _, _, rho_cool_2 = get_gas_properties(30.0, P_outlet, o2_frac, n2_frac)
        vol_flow_cool_2 = (mass_flow_2 / rho_cool_2) * SEC_PER_HOUR
        st.metric("Vol Flow after cooling to 30°C", f"{vol_flow_cool_2:.1f} m³/h")
        
        with st.expander("Stage Details - Section 2"):
            for stage in stages_data_2:
                st.write(f"**Stage {stage['stage']}**")
                st.write(f"  Pressure: {stage['P_in']:.2f} → {stage['P_out']:.2f} bar")
                st.write(f"  Temperature: {stage['T_in']:.1f} → {stage['T_out']:.1f} °C")
                st.write(f"  Gas Power: {stage['power_gas']:.3f} MW")
                st.write(f"  Shaft Power: {stage['power_shaft']:.3f} MW")
                st.write("---")
    
    # Summary
    st.header("Overall Summary")
    col3, col4, col5, col6 = st.columns(4)
    
    with col3:
        st.metric("Total Gas Power", f"{power_gas_1 + power_gas_2:.3f} MW",
                  help="Thermodynamic power transferred to gas")
    
    with col4:
        st.metric("Total Shaft Power", f"{power_shaft_1 + power_shaft_2:.3f} MW",
                  help="Total power required at motor shaft")
    
    with col5:
        st.metric("Overall Pressure Ratio", f"{P_outlet / P_inlet:.1f}")
    
    with col6:
        MW_mix, _, _, _ = get_gas_properties(T_inlet, P_inlet, o2_frac, n2_frac)
        st.metric("Gas Molecular Weight", f"{MW_mix:.2f} g/mol")
    
    # Efficiency loss breakdown
    with st.expander("Power Loss Analysis"):
        mechanical_loss_1 = power_shaft_1 - power_gas_1
        mechanical_loss_2 = power_shaft_2 - power_gas_2
        total_mechanical_loss = mechanical_loss_1 + mechanical_loss_2
        
        st.write("**Mechanical Losses:**")
        st.write(f"Section 1: {mechanical_loss_1:.3f} MW ({mechanical_loss_1/power_shaft_1*100:.1f}%)")
        st.write(f"Section 2: {mechanical_loss_2:.3f} MW ({mechanical_loss_2/power_shaft_2*100:.1f}%)")
        st.write(f"Total: {total_mechanical_loss:.3f} MW")

# Function to generate table data
def generate_table(final_pressure, eta_isen, eta_mech):
    """Generate comprehensive table for given final pressure"""
    
    inlet_flows = [100000, 200000, 300000, 400000, 600000]
    T_inlet_table = 20.0
    P_inlet_table = 1.0
    P_intermediate_table = 7.0
    stages_1_table = 2
    stages_2_table = 3
    T_cool = 30.0
    
    results = []
    
    for vol_in in inlet_flows:
        # Section 1: Compress to 7 bar
        T_out_1, power_gas_1, power_shaft_1, mass_flow_1, _ = compress_multi_stage(
            T_inlet_table, P_inlet_table, P_intermediate_table, vol_in,
            stages_1_table, o2_frac, n2_frac, eta_isen, eta_mech
        )
        
        # Volume at outlet of section 1 (before cooling)
        _, _, _, rho_out_1 = get_gas_properties(T_out_1, P_intermediate_table, o2_frac, n2_frac)
        vol_out_1 = (mass_flow_1 / rho_out_1) * SEC_PER_HOUR
        
        # Volume after cooling to 30°C
        _, _, _, rho_cool_1 = get_gas_properties(T_cool, P_intermediate_table, o2_frac, n2_frac)
        vol_cool_1 = (mass_flow_1 / rho_cool_1) * SEC_PER_HOUR
        
        # Section 2: Compress to final pressure
        T_out_2, power_gas_2, power_shaft_2, mass_flow_2, _ = compress_multi_stage(
            T_cool, P_intermediate_table, final_pressure, vol_cool_1,
            stages_2_table, o2_frac, n2_frac, eta_isen, eta_mech
        )
        
        # Volume at outlet of section 2 (before cooling)
        _, _, _, rho_out_2 = get_gas_properties(T_out_2, final_pressure, o2_frac, n2_frac)
        vol_out_2 = (mass_flow_2 / rho_out_2) * SEC_PER_HOUR
        
        # Volume after final cooling to 30°C
        _, _, _, rho_cool_2 = get_gas_properties(T_cool, final_pressure, o2_frac, n2_frac)
        vol_cool_2 = (mass_flow_2 / rho_cool_2) * SEC_PER_HOUR
        
        # Total powers
        total_power_gas = power_gas_1 + power_gas_2
        total_power_shaft = power_shaft_1 + power_shaft_2
        
        results.append({
            'Inlet Vol (m³/h)': vol_in,
            'Inlet T (°C)': T_inlet_table,
            'Inlet P (bar)': P_inlet_table,
            'Sec1 Out P (bar)': P_intermediate_table,
            'Sec1 Stages': stages_1_table,
            'Sec1 Out T (°C)': f"{T_out_1:.1f}",
            'Sec1 Out Vol (m³/h)': f"{vol_out_1:.1f}",
            'Sec1 Vol@30°C (m³/h)': f"{vol_cool_1:.1f}",
            'Sec1 Gas Power (MW)': f"{power_gas_1:.3f}",
            'Sec1 Shaft Power (MW)': f"{power_shaft_1:.3f}",
            'Sec2 Out P (bar)': final_pressure,
            'Sec2 Stages': stages_2_table,
            'Sec2 Out T (°C)': f"{T_out_2:.1f}",
            'Sec2 Out Vol (m³/h)': f"{vol_out_2:.1f}",
            'Sec2 Vol@30°C (m³/h)': f"{vol_cool_2:.1f}",
            'Sec2 Gas Power (MW)': f"{power_gas_2:.3f}",
            'Sec2 Shaft Power (MW)': f"{power_shaft_2:.3f}",
            'Total Gas Power (MW)': f"{total_power_gas:.3f}",
            'Total Shaft Power (MW)': f"{total_power_shaft:.3f}"
        })
    
    return pd.DataFrame(results)

# Generate tables for each pressure
with tab2:
    st.header("Compression Table - Final Pressure: 45 bar")
    st.info(f"Using η_isentropic = {eta_isentropic:.0%}, η_mechanical = {eta_mech:.0%}")
    df_45 = generate_table(45, eta_isentropic, eta_mech)
    st.dataframe(df_45, use_container_width=True)
    st.download_button("Download CSV", df_45.to_csv(index=False), 
                       "compression_45bar.csv", "text/csv", key="dl_45")

with tab3:
    st.header("Compression Table - Final Pressure: 50 bar")
    st.info(f"Using η_isentropic = {eta_isentropic:.0%}, η_mechanical = {eta_mech:.0%}")
    df_50 = generate_table(50, eta_isentropic, eta_mech)
    st.dataframe(df_50, use_container_width=True)
    st.download_button("Download CSV", df_50.to_csv(index=False), 
                       "compression_50bar.csv", "text/csv", key="dl_50")

with tab4:
    st.header("Compression Table - Final Pressure: 70 bar")
    st.info(f"Using η_isentropic = {eta_isentropic:.0%}, η_mechanical = {eta_mech:.0%}")
    df_70 = generate_table(70, eta_isentropic, eta_mech)
    st.dataframe(df_70, use_container_width=True)
    st.download_button("Download CSV", df_70.to_csv(index=False), 
                       "compression_70bar.csv", "text/csv", key="dl_70")

with tab5:
    st.header("Compression Table - Final Pressure: 84 bar")
    st.info(f"Using η_isentropic = {eta_isentropic:.0%}, η_mechanical = {eta_mech:.0%}")
    df_84 = generate_table(84, eta_isentropic, eta_mech)
    st.dataframe(df_84, use_container_width=True)
    st.download_button("Download CSV", df_84.to_csv(index=False), 
                       "compression_84bar.csv", "text/csv", key="dl_84")
