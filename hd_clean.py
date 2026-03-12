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