import streamlit as st
import sqlite3
import pandas as pd

st.set_page_config(page_title="Huntington Predictor", layout="wide")

# ================= DATABASE =================

conn = sqlite3.connect("patients.db", check_same_thread=False)
c = conn.cursor()

c.execute("""
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

    user = st.text_input("Username")
    pw = st.text_input("Password", type="password")

    if st.button("Login"):
        if user == "doctor" and pw == "123":
            st.session_state.login = True
        else:
            st.error("Invalid")

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

# ================= UI =================

st.title("Huntington Disease Risk Predictor")

name = st.text_input("Patient Name")
age = st.number_input("Age", 1, 100)
cag = st.slider("CAG Repeat Count", 10, 60)

if st.button("Predict Risk"):

    risk = risk_predict(cag)

    if risk == "No Risk":
        st.success(risk)
    elif risk == "Risk":
        st.warning(risk)
    else:
        st.error(risk)

    # ⭐⭐⭐ FIXED INSERT QUERY ⭐⭐⭐
    c.execute(
        "INSERT INTO patients(name,age,cag,risk) VALUES (?,?,?,?)",
        (name, age, cag, risk)
    )
    conn.commit()

# ================= DATABASE VIEW =================

st.subheader("Patient History")

df = pd.read_sql_query("SELECT * FROM patients ORDER BY id DESC", conn)

st.dataframe(df)

# ================= EXPORT LAST =================

if st.button("Export Last Patient Report"):

    if len(df) > 0:
        last = df.iloc[0]

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
