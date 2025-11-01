import streamlit as st

st.title("🚗 Simulasi Gerak Lurus Beraturan (GLB) - Mobil Interaktif")

# Baca file HTML
with open("Simulasi_Gerak_Lurus_Beraturan.html", "r", encoding="utf-8") as f:
    html_content = f.read()

# Tampilkan HTML di Streamlit
st.components.v1.html(html_content, height=900, scrolling=True)
