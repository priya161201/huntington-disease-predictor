import streamlit as st
import numpy as np
import pandas as pd
import sqlite3
import joblib
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
from sklearn.metrics import roc_curve, auc, confusion_matrix
from reportlab.pdfgen import canvas
from io import BytesIO

# ------------------ PAGE CONFIG ------------------

st.set_page_config(page_title="Genomic Risk Platform", layout="wide")

# ------------------ DATABASE ------------------

conn = sqlite3.connect("patients.db", check_same_thread=False)
c = conn.cursor()

c.execute("""
CREATE TABLE IF NOT EXISTS patients(
name TEXT,
cag INTEGER,
risk TEXT
)
""")
conn.commit()

# ------------------ LOGIN ------------------

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

# ------------------ THEME TOGGLE ------------------

theme = st.sidebar.radio("Theme", ["Light", "Dark"], key="theme_toggle")

# ------------------ MODEL ------------------

model = joblib.load("huntington_ml_model.pkl")

# ------------------ TABS ------------------

tab1, tab2, tab3, tab4, tab5 = st.tabs(
["Prediction","Database","Genome Track","Analytics","Report"]
)

# ==========================================================
# TAB 1 PREDICTION
# ==========================================================

with tab1:

    
# ================= RISK FUNCTION =================

def risk_predict(cag):
    if cag < 20:
        return "No Risk"
    elif cag <= 37:
        return "Risk"
    else:
        return "High Risk"

# ================= UI =================

st.title("Huntington Disease Risk Predictor")

name = st.text_input("Patient Name")
age = st.number_input("Age", 1, 100)
cag = st.slider("CAG Repeat Count", 10, 60)

if st.button("Predict Risk"):

    risk = risk_predict(cag)

    if risk == "No Risk":
        st.success(risk)
    elif risk == "Risk":
        st.warning(risk)
    else:
        st.error(risk)

    # ⭐⭐⭐ FIXED INSERT QUERY ⭐⭐⭐
    c.execute(
        "INSERT INTO patients(name,age,cag,risk) VALUES (?,?,?,?)",
        (name, age, cag, risk)
    )
    conn.commit()


        col1,col2 = st.columns(2)

        with col1:
            st.metric("Predicted Risk",risk)

        with col2:

            fig = go.Figure(go.Indicator(
                mode="gauge+number",
                value=prob*100,
                title={"text":"Risk %"},
                gauge={
                    "axis":{"range":[0,100]},
                    "bar":{"color":"red"},
                    "steps":[
                        {"range":[0,36],"color":"green"},
                        {"range":[36,100],"color":"red"}
                    ]
                }
            ))

            st.plotly_chart(fig,use_container_width=True)

# ==========================================================
# TAB 2 DATABASE
# ==========================================================

with tab2:

    st.header("Patient Database")

    pname = st.text_input("Patient Name",key="db_name")
    pcag = st.number_input("CAG",10,60,20,key="db_cag")

    if st.button("Save Patient",key="save_db"):

        pred = model.predict(pd.DataFrame({"CAG_Repeats":[pcag]}))[0]
        risk = "High" if pred==1 else "Low"

        c.execute("INSERT INTO patients VALUES(?,?,?)",(pname,pcag,risk))
        conn.commit()
        st.success("Saved")

    search = st.text_input("Search Patient",key="search_db")

    df = pd.read_sql("SELECT * FROM patients",conn)

    if search:
        df = df[df["name"].str.contains(search)]

    st.dataframe(df)

# ==========================================================
# TAB 3 GENOME TRACK
# ==========================================================

with tab3:

    st.header("Genome Browser Visualization")

    positions = np.random.randint(1,10000,20)
    lengths = np.random.randint(1,10,20)

    plt.figure(figsize=(10,3))
    for p,l in zip(positions,lengths):
        plt.plot([p,p+l],[1,1],linewidth=l)
    plt.yticks([])
    plt.xlabel("Genome Position")
    st.pyplot(plt)

# ==========================================================
# TAB 4 ANALYTICS
# ==========================================================

with tab4:

    st.header("Model Analytics")

    X = np.random.randint(15,55,200)
    y = np.array([1 if i>36 else 0 for i in X])

    y_prob = model.predict_proba(pd.DataFrame({"CAG_Repeats":X}))[:,1]

    fpr,tpr,_ = roc_curve(y,y_prob)
    roc_auc = auc(fpr,tpr)

    fig = plt.figure()
    plt.plot(fpr,tpr,label=f"AUC={roc_auc:.2f}")
    plt.plot([0,1],[0,1],"--")
    plt.legend()
    st.pyplot(fig)

    cm = confusion_matrix(y,model.predict(pd.DataFrame({"CAG_Repeats":X})))

    fig2 = plt.figure()
    sns.heatmap(cm,annot=True,fmt="d")
    st.pyplot(fig2)

# ==========================================================
# TAB 5 PDF REPORT
# ==========================================================

with tab5:

    st.header("Generate Medical Report")

    rname = st.text_input("Report Patient Name",key="report_name")
    rcag = st.number_input("Report CAG",10,60,20,key="report_cag")

    if st.button("Generate PDF",key="pdf_btn"):

        buffer = BytesIO()
        cpdf = canvas.Canvas(buffer)

        pred = model.predict(pd.DataFrame({"CAG_Repeats":[rcag]}))[0]
        risk = "High Risk" if pred==1 else "Low Risk"

        cpdf.drawString(100,750,f"Patient: {rname}")
        cpdf.drawString(100,720,f"CAG: {rcag}")
        cpdf.drawString(100,690,f"Risk: {risk}")

        cpdf.save()

        st.download_button(
            "Download Report",
            buffer.getvalue(),
            file_name="report.pdf"
        )
