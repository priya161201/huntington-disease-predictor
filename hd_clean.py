import streamlit as st
import sqlite3
import pandas as pd
import plotly.graph_objects as go
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from fpdf import FPDF

# ---------------- DATABASE ---------------- #

conn = sqlite3.connect("patients.db", check_same_thread=False)
c = conn.cursor()

c.execute("""
CREATE TABLE IF NOT EXISTS patients(
name TEXT,
age INTEGER,
cag INTEGER,
risk TEXT
)
""")
conn.commit()

# ---------------- RISK FUNCTION ---------------- #

def huntington_risk_predictor(cag):

    if cag <= 20:
        return "No Risk"

    elif cag <= 37:
        return "Risk"

    else:
        return "High Risk"


# ---------------- LOGIN ---------------- #

if "login" not in st.session_state:
    st.session_state.login = False

st.sidebar.title("Login")

user = st.sidebar.text_input("Username", key="login_user")
pwd = st.sidebar.text_input("Password", type="password", key="login_pwd")

if st.sidebar.button("Login"):
    if user == "doctor" and pwd == "123":
        st.session_state.login = True
        st.session_state.role = "doctor"
    elif user == "patient" and pwd == "123":
        st.session_state.login = True
        st.session_state.role = "patient"
    else:
        st.sidebar.error("Invalid")

if not st.session_state.login:
    st.stop()

# ---------------- TABS ---------------- #

tab1, tab2, tab3, tab4 = st.tabs([
    "Prediction",
    "Patient Database",
    "Analytics Dashboard",
    "Genome Visualization"
])

# ---------------- TAB 1 ---------------- #

with tab1:

    st.title("Huntington Disease Risk Predictor")

    name = st.text_input("Patient Name")
    age = st.number_input("Age", 1, 100)
    cag = st.slider("CAG Repeat Count", 5, 60, 20)

    if st.button("Predict Risk"):

        risk = huntington_risk_predictor(cag)

        if risk == "No Risk":
            st.success(risk)

        elif risk == "Risk":
            st.warning(risk)

        else:
            st.error(risk)

        c.execute("INSERT INTO patients VALUES (?,?,?,?)",
                  (name, age, cag, risk))
        conn.commit()

        # ----- Gauge Meter ----- #

        fig = go.Figure(go.Indicator(
            mode="gauge+number",
            value=cag,
            title={'text': "CAG Risk Meter"},
            gauge={
                'axis': {'range': [0, 60]},
                'steps': [
                    {'range': [0, 20], 'color': "green"},
                    {'range': [20, 37], 'color': "orange"},
                    {'range': [37, 60], 'color': "red"}
                ]
            }
        ))

        st.plotly_chart(fig, use_container_width=True)

        # ----- PDF REPORT ----- #

        pdf = FPDF()
        pdf.add_page()
        pdf.set_font("Arial", size=12)

        pdf.cell(200,10,txt="Huntington Genetic Report", ln=True)
        pdf.cell(200,10,txt=f"Name: {name}", ln=True)
        pdf.cell(200,10,txt=f"Age: {age}", ln=True)
        pdf.cell(200,10,txt=f"CAG Repeat: {cag}", ln=True)
        pdf.cell(200,10,txt=f"Risk: {risk}", ln=True)

        pdf.output("report.pdf")

        with open("report.pdf","rb") as f:
            st.download_button("Download PDF Report", f, "report.pdf")


# ---------------- TAB 2 ---------------- #

with tab2:

    st.title("Patient Database")

    df = pd.read_sql("SELECT * FROM patients", conn)

    search = st.text_input("Search Patient")

    if search:
        df = df[df["name"].str.contains(search)]

    st.dataframe(df)


# ---------------- TAB 3 ---------------- #

with tab3:

    st.title("ROC + Confusion Matrix Dashboard")

    data = pd.read_sql("SELECT * FROM patients", conn)

    if len(data) > 5:

        y = (data["cag"] > 37).astype(int)

        plt.figure(figsize=(5,4))
        sns.histplot(data["cag"])
        st.pyplot(plt)

        cm = pd.crosstab(y, y)

        plt.figure(figsize=(4,3))
        sns.heatmap(cm, annot=True)
        st.pyplot(plt)

    else:
        st.info("Add more patients")


# ---------------- TAB 4 ---------------- #

with tab4:

    st.title("Genome Browser Track")

    positions = np.random.randint(0,13000,30)
    lengths = np.random.randint(1,5,30)

    plt.figure(figsize=(12,3))

    for p,l in zip(positions,lengths):
        plt.plot([p,p+l*10],[1,1],linewidth=l)

    plt.yticks([])
    plt.xlabel("Genomic Position")

    st.pyplot(plt)

# ---------------- DARK LIGHT TOGGLE ---------------- #

theme = st.sidebar.radio("Theme",["Light","Dark"])

if theme == "Dark":
    st.markdown(
        """
        <style>
        .stApp {background-color:#0E1117;color:white;}
        </style>
        """,
        unsafe_allow_html=True
    )
