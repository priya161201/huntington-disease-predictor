import streamlit as st
import pandas as pd
import numpy as np
import joblib
import plotly.graph_objects as go
from streamlit_lottie import st_lottie
import json

# ---------------- PAGE CONFIG ----------------
st.set_page_config(page_title="Genomic AI Platform", layout="wide")

# ---------------- THEME TOGGLE ----------------
theme = st.sidebar.toggle("🌙 Dark Mode")

if theme:
    st.markdown("""
    <style>
    .stApp {background-color:#0E1117;color:white;}
    </style>
    """, unsafe_allow_html=True)

# ---------------- LOGIN SYSTEM ----------------
if "login" not in st.session_state:
    st.session_state.login = False

def login():
    st.title("🔐 Genomic AI Login")

    user = st.text_input("Username")
    pwd = st.text_input("Password", type="password")

    if st.button("Login"):
        if user == "admin" and pwd == "genomics":
            st.session_state.login = True
            st.success("Login Success")
            st.rerun()
        else:
            st.error("Invalid Credentials")

if not st.session_state.login:
    login()
    st.stop()

# ---------------- LOAD MODEL ----------------
model = joblib.load("huntington_ml_model.pkl")

# ---------------- SIDEBAR ----------------
st.sidebar.title("🧬 Navigation")

page = st.sidebar.radio(
    "Go to",
    ["🏠 Dashboard", "👤 Patient Prediction", "📊 Batch Prediction"]
)

# ---------------- ANIMATION ----------------
def load_lottie():
    return {
        "v":"5.5.7",
        "fr":30,
        "ip":0,
        "op":60,
        "w":200,
        "h":200,
        "nm":"DNA",
        "ddd":0,
        "assets":[],
        "layers":[]
    }

st_lottie(load_lottie(), height=150)

# ---------------- DATABASE ----------------
if "patients" not in st.session_state:
    st.session_state.patients = pd.DataFrame(columns=["Name","CAG","Risk"])

# ================= DASHBOARD =================
if page == "🏠 Dashboard":

    st.title("🧬 Huntington Genomic Intelligence Platform")

    col1,col2,col3 = st.columns(3)

    col1.metric("Total Patients", len(st.session_state.patients))
    col2.metric("High Risk", (st.session_state.patients["Risk"]=="High").sum())
    col3.metric("Normal", (st.session_state.patients["Risk"]=="Low").sum())

    # SEARCH FILTER
    search = st.text_input("🔍 Search Patient")

    if search:
        df = st.session_state.patients[
            st.session_state.patients["Name"].str.contains(search, case=False)
        ]
        st.dataframe(df)

    else:
        st.dataframe(st.session_state.patients)

# ================= SINGLE PREDICTION =================
if page == "👤 Patient Prediction":

    st.title("🧠 Real-Time Huntington Risk Prediction")

    name = st.text_input("Patient Name")

    cag = st.slider("CAG Repeat Count", 10, 60, 25)

    if st.button("Predict"):

        pred = model.predict(pd.DataFrame({"CAG_Repeats":[cag]}))[0]

        risk = "High" if pred==1 else "Low"

        st.success(f"Risk Level: {risk}")

        # GAUGE METER
        fig = go.Figure(go.Indicator(
            mode="gauge+number",
            value=cag,
            title={'text': "CAG Risk Meter"},
            gauge={
                'axis': {'range': [10,60]},
                'steps':[
                    {'range':[10,26],'color':'green'},
                    {'range':[27,35],'color':'orange'},
                    {'range':[36,60],'color':'red'}
                ]
            }
        ))

        st.plotly_chart(fig, use_container_width=True)

        # SAVE PATIENT
        new = pd.DataFrame({
            "Name":[name],
            "CAG":[cag],
            "Risk":[risk]
        })

        st.session_state.patients = pd.concat(
            [st.session_state.patients,new],
            ignore_index=True
        )

# ================= BATCH PREDICTION =================
if page == "📊 Batch Prediction":

    st.title("📁 Upload Patient CSV")

    file = st.file_uploader("Upload CSV")

    if file:

        df = pd.read_csv(file)

        df["Risk"] = model.predict(df)

        st.dataframe(df)

        st.download_button(
            "Download Results",
            df.to_csv(index=False),
            "prediction.csv"
        )
