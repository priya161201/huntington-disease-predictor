import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from io import StringIO

# ================= THEME =================

st.set_page_config(page_title="Huntington AI", layout="wide")

if "theme" not in st.session_state:
    st.session_state.theme = "Light"

if st.sidebar.radio("Theme", ["Light","Dark"]) == "Dark":
    st.session_state.theme = "Dark"
    st.markdown(
        """
        <style>
        .stApp {background-color:#0E1117;color:white}
        </style>
        """,
        unsafe_allow_html=True
    )

# ================= LOGIN =================

if "logged_in" not in st.session_state:
    st.session_state.logged_in = False

st.sidebar.title("Login")

user = st.sidebar.text_input("Username")
pw = st.sidebar.text_input("Password", type="password")

if st.sidebar.button("Login"):
    if user != "" and pw != "":
        st.session_state.logged_in = True
        st.success("Login Successful")

# ================= MAIN =================

if st.session_state.logged_in:

    st.title("🧬 Huntington Disease AI Platform")

    tab1, tab2, tab3 = st.tabs([
        "Prediction",
        "Batch Patients",
        "FASTA Analysis"
    ])

    # ================= PREDICTION =================

    with tab1:

        cag = st.slider("Select CAG Repeat", 10, 60, 20)

        if st.button("Predict Risk"):

            if cag < 27:
                risk = 10
                label = "Normal"
            elif cag < 36:
                risk = 50
                label = "Intermediate"
            else:
                risk = 90
                label = "High Risk"

            fig = go.Figure(go.Indicator(
                mode="gauge+number",
                value=risk,
                title={'text': label},
                gauge={
                    'axis': {'range': [0,100]},
                    'bar': {'color': "red"},
                    'steps': [
                        {'range': [0,30], 'color': "green"},
                        {'range': [30,70], 'color': "orange"},
                        {'range': [70,100], 'color': "red"},
                    ]
                }
            ))

            st.plotly_chart(fig, use_container_width=True)

    # ================= BATCH =================

    with tab2:

        st.subheader("Upload Patient CSV")

        file = st.file_uploader("Upload CSV", type=["csv"])

        if file:

            df = pd.read_csv(file)

            def predict(x):
                if x < 27:
                    return "Normal"
                elif x < 36:
                    return "Intermediate"
                else:
                    return "High Risk"

            df["Prediction"] = df["CAG_Repeats"].apply(predict)

            st.dataframe(df)

            st.download_button(
                "Download Results",
                df.to_csv(index=False),
                "predictions.csv"
            )

    # ================= FASTA =================

    with tab3:

        fasta = st.file_uploader("Upload FASTA", type=["fasta","fa"])

        if fasta:

            stringio = StringIO(fasta.getvalue().decode("utf-8"))
            seq = stringio.read()

            cag_count = seq.count("CAG")

            st.metric("Detected CAG Motifs", cag_count)

else:

    st.warning("Please login from sidebar")
