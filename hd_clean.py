import streamlit as st
import joblib
import pandas as pd
import matplotlib.pyplot as plt

st.set_page_config(page_title="Huntington Predictor", layout="centered")

model = joblib.load("huntington_ml_model.pkl")

tab1, tab2 = st.tabs(["🧬 Prediction", "📚 About Huntington Disease"])

with tab1:

    st.title("Huntington Disease Risk Predictor")

    repeat = st.slider("Select CAG Repeat Count", 10, 60, 20)

    if st.button("Predict Risk"):

        result = model.predict(
            pd.DataFrame({"CAG_Repeats":[repeat]})
        )[0]

        if repeat < 27:
            st.success("Normal Range")
        elif repeat < 36:
            st.warning("Intermediate Risk")
        else:
            st.error("Huntington Disease Risk")

        fig, ax = plt.subplots()
        ax.bar(["Your Repeat"], [repeat])
        ax.axhline(36, color="red", linestyle="--")
        ax.set_ylabel("CAG Count")

        st.pyplot(fig)

with tab2:

    st.header("What is Huntington Disease?")

    st.write("""
    Huntington’s Disease is a genetic neurodegenerative disorder caused by 
    expansion of CAG trinucleotide repeats in the HTT gene.

    • Normal: < 27 repeats  
    • Intermediate: 27–35 repeats  
    • Disease Risk: ≥ 36 repeats  

    This tool demonstrates computational genomics + AI based risk prediction.

    """)

from Bio import SeqIO
import re

tab1, tab2, tab3 = st.tabs([
    "🧬 Prediction",
    "📚 About Disease",
    "🔬 FASTA Analysis"
])

with tab3:

    st.header("Upload HTT Gene FASTA")

    fasta_file = st.file_uploader("Upload FASTA file", type=["fasta","fa"])

    if fasta_file:

        record = SeqIO.read(fasta_file, "fasta")
        seq = str(record.seq)

        repeats = re.findall(r"(?:CAG){2,}", seq)
        counts = [len(r)//3 for r in repeats]

        if counts:
            max_repeat = max(counts)
            st.write("Longest CAG Repeat:", max_repeat)

            if max_repeat >= 36:
                st.error("⚠️ Huntington Risk Detected")
            else:
                st.success("Normal Repeat Range")
        else:
            st.info("No significant repeat region found")

from sklearn.metrics import roc_curve, auc
import numpy as np

with tab4:

    st.header("Model ROC Curve")

    x = np.linspace(15,55,100)
    y = model.predict_proba(
        pd.DataFrame({"CAG_Repeats":x})
    )[:,1]

    fpr, tpr, _ = roc_curve((x>=36).astype(int), y)
    roc_auc = auc(fpr,tpr)

    fig, ax = plt.subplots()
    ax.plot(fpr,tpr,label=f"AUC={roc_auc:.2f}")
    ax.plot([0,1],[0,1],'--')
    ax.legend()
    st.pyplot(fig)

st.header("Batch Prediction")

csv = st.file_uploader("Upload Patient CSV", type=["csv"])

if csv:

    df = pd.read_csv(csv)

    df["Prediction"] = model.predict(df[["CAG_Repeats"]])

    st.dataframe(df)

    st.download_button(
        "Download Results",
        df.to_csv(index=False),
        "predictions.csv"
    )

st.set_page_config(
    page_title="Genomic Disease Dashboard",
    page_icon="🧬",
    layout="wide"
)

st.markdown("""
<style>
.big-font {
font-size:40px !important;
font-weight:700;
color:#2E86C1;
}
</style>
""", unsafe_allow_html=True)

st.markdown('<p class="big-font">Genomic Risk Intelligence System</p>', unsafe_allow_html=True)
