import streamlit as st
import pandas as pd
import joblib
import matplotlib.pyplot as plt
import numpy as np
import re
import io
from Bio import SeqIO
from sklearn.metrics import roc_curve, auc

# ---------------- PAGE CONFIG ---------------- #

st.set_page_config(
    page_title="Genomic Risk Intelligence System",
    page_icon="🧬",
    layout="wide"
)

st.markdown("""
<style>
.big-font {
font-size:42px !important;
font-weight:700;
color:#2E86C1;
}
</style>
""", unsafe_allow_html=True)

st.markdown(
'<p class="big-font">Genomic Risk Intelligence Dashboard</p>',
unsafe_allow_html=True
)

# ---------------- LOAD MODEL ---------------- #

model = joblib.load("huntington_ml_model.pkl")

# ---------------- TABS ---------------- #

tab1, tab2, tab3, tab4 = st.tabs([
    "🧬 Single Prediction",
    "📂 Batch Prediction",
    "🔬 FASTA Genomic Analysis",
    "📊 Model ROC Dashboard"
])

# ==========================================================
# 🧬 SINGLE PREDICTION
# ==========================================================

with tab1:

    st.header("Huntington Disease Risk Predictor")

    repeat = st.slider("Select CAG Repeat Count", 10, 60, 20)

    if st.button("Predict Risk"):

        pred = model.predict(
            pd.DataFrame({"CAG_Repeats":[repeat]})
        )[0]

        if repeat < 27:
            st.success("✅ Normal Range")
        elif repeat < 36:
            st.warning("⚠️ Intermediate Risk")
        else:
            st.error("🚨 Huntington Disease Risk")

        fig, ax = plt.subplots()

        ax.bar(["Your Repeat"], [repeat], color="skyblue")
        ax.axhline(36, color="red", linestyle="--",
                   label="Disease Threshold")

        ax.set_ylabel("CAG Count")
        ax.legend()

        st.pyplot(fig)

# ==========================================================
# 📂 BATCH PREDICTION
# ==========================================================

with tab2:

    st.header("Batch Patient Prediction")

    st.info("Upload CSV containing column: CAG_Repeats")

    file = st.file_uploader("Upload CSV", type=["csv"])

    if file:

        df = pd.read_csv(file)

        df["Prediction"] = model.predict(
            df[["CAG_Repeats"]]
        )

        st.dataframe(df)

        st.download_button(
            "Download Results",
            df.to_csv(index=False),
            "huntington_predictions.csv"
        )

# ==========================================================
# 🔬 FASTA GENOMIC ANALYSIS (FIXED CLOUD VERSION)
# ==========================================================

with tab3:

    st.header("HTT Gene FASTA Analysis")

    fasta = st.file_uploader(
        "Upload FASTA file", type=["fasta","fa"]
    )

    if fasta:

        text_stream = io.StringIO(
            fasta.getvalue().decode("utf-8")
        )

        record = SeqIO.read(text_stream, "fasta")

        seq = str(record.seq)

        cag_blocks = re.findall(r"(?:CAG){2,}", seq)

        counts = [len(b)//3 for b in cag_blocks]

        if counts:

            max_repeat = max(counts)

            st.subheader(
                f"Longest CAG Repeat Found: {max_repeat}"
            )

            if max_repeat >= 36:
                st.error(
                "⚠️ Huntington Disease Risk Region Detected"
                )
            else:
                st.success("Normal Repeat Length")

        else:
            st.info(
            "No significant CAG repeat expansion found"
            )

# ==========================================================
# 📊 ROC CURVE DASHBOARD
# ==========================================================

with tab4:

    st.header("Model ROC Curve")

    x = np.linspace(15,55,200)

    probs = model.predict_proba(
        pd.DataFrame({"CAG_Repeats":x})
    )[:,1]

    y_true = (x>=36).astype(int)

    fpr, tpr, _ = roc_curve(y_true, probs)

    roc_auc = auc(fpr, tpr)

    fig, ax = plt.subplots()

    ax.plot(fpr, tpr,
            label=f"AUC = {roc_auc:.2f}",
            color="darkorange")

    ax.plot([0,1],[0,1],'--', color="navy")

    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title("ROC Curve")

    ax.legend()

    st.pyplot(fig)
