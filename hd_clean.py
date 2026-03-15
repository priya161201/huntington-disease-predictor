import streamlit as st
import numpy as np
import pandas as pd
import joblib
import sqlite3
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from Bio import SeqIO

# ---------------------- PAGE CONFIG ----------------------

st.set_page_config(
    page_title="Genomic AI Platform",
    layout="wide"
)

# ---------------------- THEME TOGGLE ----------------------

theme = st.sidebar.radio("Theme", ["Light", "Dark"])

if theme == "Dark":
    st.markdown("""
    <style>
    body {background-color:#0E1117;color:white;}
    </style>
    """, unsafe_allow_html=True)

# ---------------------- LOGIN SYSTEM ----------------------

st.sidebar.title("Login")

username = st.sidebar.text_input("Username")
password = st.sidebar.text_input("Password", type="password")

if username != "admin" or password != "1234":
    st.warning("Login to continue")
    st.stop()

# ---------------------- DATABASE ----------------------

conn = sqlite3.connect("patients.db", check_same_thread=False)
c = conn.cursor()

c.execute("""
CREATE TABLE IF NOT EXISTS patients(
id INTEGER PRIMARY KEY AUTOINCREMENT,
name TEXT,
cag INTEGER,
risk TEXT
)
""")
conn.commit()

# ---------------------- LOAD MODEL ----------------------

try:
    model = joblib.load("huntington_ml_model.pkl")
except:
    model = None

# ---------------------- TABS ----------------------

tab1, tab2, tab3, tab4 = st.tabs([
    "🧠 Risk Prediction",
    "🧬 FASTA Analysis",
    "📊 Analytics",
    "🗄 Patient Database"
])

# =========================================================
# TAB 1 — RISK PREDICTION
# =========================================================

with tab1:

    st.title("🧠 Huntington Disease Risk Predictor")

    cag = st.slider("CAG Repeat Count", 10, 60, 20)

    if st.button("Predict Risk"):

        if model is not None:
            pred = model.predict(pd.DataFrame({"CAG_Repeats":[cag]}))[0]
        else:
            pred = 1 if cag >= 36 else 0

        risk = "High Risk" if pred == 1 else "Low Risk"

        st.success(f"Risk Status: {risk}")

        # Gauge meter
        fig = go.Figure(go.Indicator(
            mode="gauge+number",
            value=cag,
            title={"text":"CAG Severity"},
            gauge={
                "axis":{"range":[10,60]},
                "steps":[
                    {"range":[10,27],"color":"green"},
                    {"range":[27,36],"color":"yellow"},
                    {"range":[36,60],"color":"red"}
                ]
            }
        ))

        st.plotly_chart(fig, use_container_width=True)

        name = st.text_input("Patient Name")

        if st.button("Save Patient"):
            c.execute(
                "INSERT INTO patients(name,cag,risk) VALUES(?,?,?)",
                (name,cag,risk)
            )
            conn.commit()
            st.success("Saved")

# =========================================================
# TAB 2 — FASTA ANALYSIS
# =========================================================

with tab2:

    st.title("🧬 HTT FASTA Analysis")

    fasta = st.file_uploader("Upload FASTA", type=["fa","fasta"])

    if fasta is not None:

        import io
        record = SeqIO.read(io.StringIO(fasta.getvalue().decode()), "fasta")

        seq = str(record.seq)

        count = seq.count("CAG")

        st.metric("Total CAG Count", count)

        lengths = [len(x) for x in seq.split("CAG")]

        plt.hist(lengths)
        st.pyplot(plt)

# =========================================================
# TAB 3 — ANALYTICS
# =========================================================

with tab3:

    st.title("📊 Population Analytics")

    data = pd.DataFrame({
        "CAG": np.random.randint(15,55,200)
    })

    data["Risk"] = data["CAG"].apply(lambda x:1 if x>=36 else 0)

    fig, ax = plt.subplots()
    ax.scatter(data["CAG"], data["Risk"])
    ax.axvline(36,color="red")
    st.pyplot(fig)

# =========================================================
# TAB 4 — DATABASE
# =========================================================

with tab4:

    st.title("🗄 Patient Records")

    df = pd.read_sql("SELECT * FROM patients", conn)

    search = st.text_input("Search Patient")

    if search:
        df = df[df["name"].str.contains(search)]

    st.dataframe(df)

    if st.button("Download CSV"):
        st.download_button(
            "Download",
            df.to_csv().encode(),
            "patients.csv"
        )
