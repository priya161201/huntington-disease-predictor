import streamlit as st
import sqlite3
from fpdf import FPDF
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt

# ---------------- DB ----------------

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

# ---------------- LOGIN ----------------

st.sidebar.title("Login")

user = st.sidebar.text_input("Username")
pw = st.sidebar.text_input("Password", type="password")

if st.sidebar.button("Login"):
    if user == "doctor" and pw == "123":
        st.session_state.login = True
    else:
        st.sidebar.error("Invalid")

if "login" not in st.session_state:
    st.stop()

# ---------------- TABS ----------------

tab1, tab2, tab3, tab4 = st.tabs([
    "Prediction",
    "Patient Database",
    "Analytics Dashboard",
    "Genome Visualization"
])

# ---------------- RISK FUNCTION ----------------

def risk_func(cag):
    if cag <= 20:
        return "No Risk"
    elif cag <= 37:
        return "Risk"
    else:
        return "High Risk"

# ---------------- TAB 1 ----------------

with tab1:

    st.title("Huntington Disease Risk Predictor")

    name = st.text_input("Patient Name")
    age = st.number_input("Age", 1, 120, 24)

    cag = st.slider("CAG Repeat Count", 10, 60, 20)

    if st.button("Predict Risk"):

        risk = risk_func(cag)

        if risk == "No Risk":
            st.success(risk)
        elif risk == "Risk":
            st.warning(risk)
        else:
            st.error(risk)

        # ⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐
        # ONLY FIXED INSERT QUERY
        # ⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐⭐

        c.execute(
            "INSERT INTO patients(name,age,cag,risk) VALUES (?,?,?,?)",
            (name, age, cag, risk)
        )
        conn.commit()

# ---------------- TAB 2 ----------------

with tab2:

    st.header("Patient History")

    df = pd.read_sql("SELECT * FROM patients", conn)

    search = st.text_input("Search Patient")

    if search:
        df = df[df["name"].str.contains(search, case=False)]

    st.dataframe(df)

# ---------------- TAB 3 ----------------

with tab3:

    st.header("Analytics Dashboard")

    df = pd.read_sql("SELECT * FROM patients", conn)

    if len(df) > 0:

        fig = go.Figure(go.Indicator(
            mode="gauge+number",
            value=df["cag"].mean(),
            title={'text': "Average CAG"},
            gauge={'axis': {'range': [10, 60]}}
        ))

        st.plotly_chart(fig)

        fig2, ax = plt.subplots()
        sns.histplot(df["cag"], bins=10, ax=ax)
        st.pyplot(fig2)

# ---------------- TAB 4 ----------------

with tab4:

    st.header("Genome Track")

    positions = np.random.randint(100, 10000, 30)
    lengths = np.random.randint(2, 25, 30)

    fig, ax = plt.subplots(figsize=(10,2))

    for p,l in zip(positions,lengths):
        ax.plot([p,p+l],[1,1], linewidth=l)

    st.pyplot(fig)

# ---------------- PDF ----------------

if st.button("Export Last Patient Report"):

    df = pd.read_sql("SELECT * FROM patients ORDER BY id DESC LIMIT 1", conn)

    if len(df) > 0:

        row = df.iloc[0]

        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("Arial", size=14)

        pdf.cell(200,10,"Huntington Genomic Report", ln=True)
        pdf.cell(200,10,f"Name: {row['name']}", ln=True)
        pdf.cell(200,10,f"Age: {row['age']}", ln=True)
        pdf.cell(200,10,f"CAG: {row['cag']}", ln=True)
        pdf.cell(200,10,f"Risk: {row['risk']}", ln=True)

        pdf.output("report.pdf")

        with open("report.pdf","rb") as f:
            st.download_button("Download PDF", f, file_name="report.pdf")
