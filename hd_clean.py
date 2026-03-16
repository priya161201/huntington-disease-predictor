import streamlit as st
import sqlite3
import pandas as pd

st.set_page_config(page_title="Huntington Predictor", layout="wide")

# ================= DATABASE =================

conn = sqlite3.connect("patients.db", check_same_thread=False)
cursor = conn.cursor()

cursor.execute("""
CREATE TABLE IF NOT EXISTS patients(
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT,
    age INTEGER,
    cag INTEGER,
    risk TEXT
)
""")

conn.commit()

# ================= LOGIN =================

if "login" not in st.session_state:
    st.session_state.login = False

with st.sidebar:
    st.title("Login")

    username = st.text_input("Username")
    password = st.text_input("Password", type="password")

    if st.button("Login"):
        if username == "doctor" and password == "123":
            st.session_state.login = True
        else:
            st.error("Invalid Login")

if not st.session_state.login:
    st.stop()

# ================= RISK FUNCTION =================

def risk_predict(cag):
    if cag < 20:
        return "No Risk"
    elif cag <= 37:
        return "Risk"
    else:
        return "High Risk"

# ================= MAIN UI =================

st.title("Huntington Disease Risk Predictor")

name = st.text_input("Patient Name")
age = st.number_input("Age", min_value=1, max_value=100, value=25)
cag = st.slider("CAG Repeat Count", 10, 60, 20)

if st.button("Predict Risk"):

    risk = risk_predict(cag)

    if risk == "No Risk":
        st.success(risk)
    elif risk == "Risk":
        st.warning(risk)
    else:
        st.error(risk)

    # INSERT INTO DATABASE
    cursor.execute(
        "INSERT INTO patients(name, age, cag, risk) VALUES (?, ?, ?, ?)",
        (name, age, cag, risk)
    )
    conn.commit()

    st.success("Patient Saved")

# ================= DATABASE TABLE =================

st.subheader("Patient History")

data = pd.read_sql_query(
    "SELECT * FROM patients ORDER BY id DESC",
    conn
)

st.dataframe(data, use_container_width=True)

# ================= EXPORT =================

if st.button("Export Last Patient Report"):

    if len(data) > 0:
        last = data.iloc[0]

        report = f"""
Patient Name : {last['name']}
Age : {last['age']}
CAG Repeat : {last['cag']}
Risk : {last['risk']}
"""

        st.download_button(
            "Download Report",
            report,
            file_name="patient_report.txt"
        )
    else:
        st.warning("No data yet")
