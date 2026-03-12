import streamlit as st
import pandas as pd
import joblib
import matplotlib.pyplot as plt
import numpy as np
import re
import io
import sqlite3
from Bio import SeqIO
from sklearn.metrics import roc_curve, auc
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

# ---------------- PAGE CONFIG ---------------- #

st.set_page_config(
    page_title="Genomic Clinical Intelligence Platform",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 Genomic Risk Intelligence Dashboard")

# ---------------- LOAD MODEL ---------------- #

model = joblib.load("huntington_ml_model.pkl")

# ---------------- DATABASE ---------------- #

conn = sqlite3.connect("patients.db")
cursor = conn.cursor()

cursor.execute("""
CREATE TABLE IF NOT EXISTS patients
(name TEXT, gene TEXT, repeat INT, result TEXT)
""")
conn.commit()

# ---------------- TABS ---------------- #

tab1, tab2, tab3, tab4 = st.tabs([
    "🧬 Single Prediction",
    "📂 Batch Prediction",
    "🔬 FASTA Genomic Analysis",
    "📊 ROC Dashboard"
])

# ==========================================================
# 🧬 SINGLE PREDICTION
# ==========================================================

with tab1:

    gene = st.selectbox(
        "Select Gene",
        ["HTT","BRCA1","CFTR"]
    )

    repeat = st.slider("CAG / Mutation Repeat", 5, 60, 20)

    name = st.text_input("Patient Name")

    if gene == "HTT":
        threshold = 36
    elif gene == "BRCA1":
        threshold = 10
    else:
        threshold = 15

    if st.button("Predict Risk"):

        if repeat >= threshold:
            result = "Risk"
            st.error("🚨 Disease Risk Detected")
        else:
            result = "Normal"
            st.success("✅ Normal Range")

        fig, ax = plt.subplots()
        ax.bar(["Repeat"], [repeat])
        ax.axhline(threshold, color="red", linestyle="--")
        st.pyplot(fig)

        # SAVE DATABASE
        cursor.execute(
            "INSERT INTO patients VALUES (?,?,?,?)",
            (name, gene, repeat, result)
        )
        conn.commit()

        st.success("Patient Saved to Database")

        # PDF REPORT
        file = "mutation_report.pdf"
        c = canvas.Canvas(file, pagesize=letter)
        c.drawString(100,750,"Genomic Mutation Report")
        c.drawString(100,700,f"Patient: {name}")
        c.drawString(100,650,f"Gene: {gene}")
        c.drawString(100,600,f"Repeat: {repeat}")
        c.drawString(100,550,f"Result: {result}")
        c.save()

        with open(file,"rb") as f:
            st.download_button(
                "Download PDF Report",
                f,
                file_name="mutation_report.pdf"
            )

    if st.button("View Patients"):

        df = pd.read_sql(
            "SELECT * FROM patients",
            conn
        )

        st.dataframe(df)

# ==========================================================
# 📂 BATCH PREDICTION
# ==========================================================

with tab2:

    st.header("Batch Prediction")

    csv = st.file_uploader("Upload CSV", type=["csv"])

    if csv:

        df = pd.read_csv(csv)

        df["Prediction"] = (
            df["CAG_Repeats"] >= 36
        ).astype(int)

        st.dataframe(df)

        st.download_button(
            "Download Results",
            df.to_csv(index=False),
            "batch_results.csv"
        )

# ==========================================================
# 🔬 FASTA ANALYSIS + VISUALIZATION
# ==========================================================

with tab3:

    fasta = st.file_uploader(
        "Upload FASTA",
        type=["fasta","fa"]
    )

    if fasta:

        text_stream = io.StringIO(
            fasta.getvalue().decode("utf-8")
        )

        record = SeqIO.read(text_stream,"fasta")

        seq = str(record.seq)

        matches = list(
            re.finditer(r"(?:CAG){2,}", seq)
        )

        positions = [
            (m.start(), (m.end()-m.start())//3)
            for m in matches
        ]

        if positions:

            pos_df = pd.DataFrame(
                positions,
                columns=["Position","Repeat_Count"]
            )

            st.dataframe(pos_df)

            fig, ax = plt.subplots()
            ax.scatter(
                pos_df["Position"],
                pos_df["Repeat_Count"],
                color="red"
            )
            ax.set_xlabel("Genome Position")
            ax.set_ylabel("Repeat Count")

            st.pyplot(fig)

        else:
            st.info("No Repeat Expansion Found")

# ==========================================================
# 📊 ROC DASHBOARD
# ==========================================================

with tab4:

    x = np.linspace(15,55,200)

    probs = model.predict_proba(
        pd.DataFrame({"CAG_Repeats":x})
    )[:,1]

    y_true = (x>=36).astype(int)

    fpr, tpr, _ = roc_curve(y_true, probs)

    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots()
    ax.plot(fpr,tpr,label=f"AUC={roc_auc:.2f}")
    ax.plot([0,1],[0,1],'--')
    ax.legend()

    st.pyplot(fig)
