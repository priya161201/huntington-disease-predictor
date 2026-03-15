import streamlit as st
import pandas as pd
import numpy as np
import sqlite3
import plotly.graph_objects as go
import plotly.express as px
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, confusion_matrix
from lifelines import KaplanMeierFitter
from reportlab.pdfgen import canvas

st.set_page_config(page_title="Genomic AI Platform", layout="wide")

# ================= DATABASE =================

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

# ================= LOGIN =================

st.sidebar.title("Login")

role = st.sidebar.selectbox("Select Role", ["Patient", "Doctor"])
username = st.sidebar.text_input("Username")
password = st.sidebar.text_input("Password", type="password")

if st.sidebar.button("Login"):
    st.session_state.logged = True
    st.session_state.role = role

if "logged" not in st.session_state:
    st.warning("Please Login from Sidebar")
    st.stop()

# ================= TABS =================

tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
"Risk Prediction",
"Patient Database",
"ROC Dashboard",
"Genome Browser",
"Survival Analysis",
"PDF Report"
])

# ================= TAB 1 =================

with tab1:

    st.title("Huntington Disease Risk Gauge")

    name = st.text_input("Patient Name")
    age = st.number_input("Age", 1, 100)
    cag = st.slider("CAG Repeat Count", 10, 60, 20)

    if st.button("Predict Risk"):

        if cag < 27:
            risk = "Normal"
            gauge_val = 20
        elif cag < 36:
            risk = "Intermediate"
            gauge_val = 50
        else:
            risk = "High Risk"
            gauge_val = 90

        fig = go.Figure(go.Indicator(
            mode="gauge+number",
            value=gauge_val,
            title={'text': risk},
            gauge={
                'axis': {'range': [0, 100]},
                'steps': [
                    {'range': [0, 30], 'color': "green"},
                    {'range': [30, 70], 'color': "orange"},
                    {'range': [70, 100], 'color': "red"},
                ]
            }
        ))

        st.plotly_chart(fig, use_container_width=True)

        if st.session_state.role == "Doctor":
            c.execute("INSERT INTO patients VALUES (?,?,?,?)",
                      (name, age, cag, risk))
            conn.commit()
            st.success("Patient Saved")

# ================= TAB 2 =================

with tab2:

    st.title("Search Patient History")

    df = pd.read_sql("SELECT * FROM patients", conn)

    search = st.text_input("Search Name")

    if search:
        df = df[df["name"].str.contains(search, case=False)]

    st.dataframe(df)

# ================= TAB 3 =================

with tab3:

    st.title("ROC + Confusion Matrix")

    X = np.random.randint(15, 55, 200).reshape(-1, 1)
    y = (X > 36).astype(int)

    model = LogisticRegression().fit(X, y)
    prob = model.predict_proba(X)[:, 1]

    fpr, tpr, _ = roc_curve(y, prob)

    roc_fig = px.line(x=fpr, y=tpr, title="ROC Curve")
    st.plotly_chart(roc_fig, use_container_width=True)

    cm = confusion_matrix(y, model.predict(X))

    cm_fig = px.imshow(cm, text_auto=True, title="Confusion Matrix")
    st.plotly_chart(cm_fig, use_container_width=True)

# ================= TAB 4 =================

with tab4:

    st.title("Genome Track Visualization")

    pos = np.random.randint(1, 10000, 30)
    length = np.random.randint(1, 20, 30)

    fig = go.Figure()

    for p, l in zip(pos, length):
        fig.add_shape(type="line",
                      x0=p, x1=p+l,
                      y0=1, y1=1,
                      line=dict(width=l))

    fig.update_layout(title="CAG Repeat Track")

    st.plotly_chart(fig, use_container_width=True)

# ================= TAB 5 =================

with tab5:

    st.title("Kaplan Survival Plot")

    T = np.random.exponential(10, 100)
    E = np.random.binomial(1, 0.6, 100)

    km = KaplanMeierFitter()
    km.fit(T, E)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=km.survival_function_.index,
        y=km.survival_function_["KM_estimate"],
        mode="lines"
    ))

    st.plotly_chart(fig, use_container_width=True)

# ================= TAB 6 =================

with tab6:

    st.title("Generate Medical PDF Report")

    pname = st.text_input("Patient Name")

    if st.button("Generate Report"):

        file = "genomic_report.pdf"

        cpdf = canvas.Canvas(file)
        cpdf.drawString(100, 750, "Genomic AI Medical Report")
        cpdf.drawString(100, 700, f"Patient: {pname}")
        cpdf.save()

        with open(file, "rb") as f:
            st.download_button("Download Report", f, file)
