import streamlit as st
import sqlite3

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

    # ✅ FIXED SQLITE INSERT
    c.execute(
        "INSERT INTO patients(name,age,cag,risk) VALUES (?,?,?,?)",
        (name, age, cag, risk)
    )

    conn.commit()

# ================= SHOW DATABASE =================

st.subheader("Patient Database")

data = c.execute("SELECT * FROM patients").fetchall()

for row in data:
    st.write(row)
