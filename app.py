import streamlit as st
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# Judul aplikasi
st.title("🧪 Simulasi Dua Fermion 1D dalam Sumur Tak Hingga")

st.markdown("""
Aplikasi ini memvisualisasikan fungsi gelombang dua fermion (spin-singlet) di sumur tak hingga 1D,
dengan interaksi delta (kontak) yang dapat diatur kekuatannya.
""")

# -------------------------
# Input parameter dari user
L = st.number_input("Panjang sumur (L)", min_value=0.5, max_value=5.0, value=1.0, step=0.1)
Nmax = st.slider("Jumlah basis (Nmax)", 2, 10, 6)
g = st.slider("Kekuatan interaksi (g)", -100.0, 100.0, 60.0)
nx = st.slider("Resolusi grid (nx)", 50, 200, 100)
# -------------------------

# Tombol untuk jalankan simulasi
if st.button("🔁 Jalankan Simulasi"):
    # Single-particle eigenfunctions dan energi (m = 1, hbar = 1 unit)
    def phi(n, x):
        return np.sqrt(2.0 / L) * np.sin(n * np.pi * x / L)

    def eps(n):
        return (np.pi**2 * n**2) / (2.0 * L**2)

    # Buat daftar pasangan (n <= m) untuk basis simetris spasial
    pairs = [(n, m) for n in range(1, Nmax+1) for m in range(n, Nmax+1)]
    B = len(pairs)

    # Precompute phi pada grid halus untuk integrasi numerik
    x_int = np.linspace(0, L, 1001)
    phi_grid = np.array([phi(n, x_int) for n in range(1, Nmax+1)])  # shape (Nmax, len(x_int))

    # Integral overlap
    def overlap(n, m, p, q):
        arr = phi_grid[n-1] * phi_grid[m-1] * phi_grid[p-1] * phi_grid[q-1]
        return np.trapz(arr, x_int)

    # Elemen matriks interaksi untuk basis simetris
    def interaction_matrix_element(pairA, pairB):
        n, m = pairA
        p, q = pairB
        I1 = overlap(n, m, p, q)
        I2 = overlap(m, n, p, q)
        I3 = overlap(n, m, q, p)
        I4 = overlap(m, n, q, p)
        return g * 0.5 * (I1 + I2 + I3 + I4)

    # Bangun Hamiltonian H = T + V pada basis simetris
    H = np.zeros((B, B))
    for i, (n, m) in enumerate(pairs):
        H[i, i] = eps(n) + eps(m)
        for j, (p, q) in enumerate(pairs):
            H[i, j] += interaction_matrix_element((n, m), (p, q))

    # Diagonalize
    Evals, Evecs = linalg.eigh(H)
    idx = np.argsort(Evals)
    Evals = Evals[idx]
    Evecs = Evecs[:, idx]

    st.success(f"Simulasi selesai! Energi keadaan dasar E₀ = {Evals[0]:.6f}")

    # Rekonstruksi fungsi gelombang keadaan dasar
    x = np.linspace(0, L, nx)
    X1, X2 = np.meshgrid(x, x, indexing='xy')
    pair_funcs = np.zeros((B, nx, nx))
    for idx_pair, (n, m) in enumerate(pairs):
        phi_n_x1 = phi(n, x)[:, None]
        phi_m_x2 = phi(m, x)[None, :]
        if n == m:
            psi_nm = phi_n_x1 @ phi_m_x2
        else:
            phi_m_x1 = phi(m, x)[:, None]
            phi_n_x2 = phi(n, x)[None, :]
            psi_nm = (phi_n_x1 @ phi_m_x2 + phi_m_x1 @ phi_n_x2) / np.sqrt(2.0)
        pair_funcs[idx_pair] = psi_nm

    c0 = Evecs[:, 0]
    Psi0 = np.tensordot(c0, pair_funcs, axes=(0, 0))
    prob = np.abs(Psi0)**2
    norm = np.trapz(np.trapz(prob, x, axis=0), x, axis=0)
    Psi0 /= np.sqrt(norm)
    prob = np.abs(Psi0)**2

    # Plot 1: kontur |Ψ|²
    fig1, ax1 = plt.subplots(figsize=(6, 5))
    cp = ax1.contourf(X1, X2, prob, levels=50)
    ax1.set_xlabel("x₁")
    ax1.set_ylabel("x₂")
    ax1.set_title(f"|Ψ(x₁,x₂)|², g={g}")
    fig1.colorbar(cp)
    st.pyplot(fig1)

    # Plot 2: densitas marginal ρ(x)
    rho_x = np.trapz(prob, x, axis=1)
    fig2, ax2 = plt.subplots(figsize=(6, 3.5))
    ax2.plot(x, rho_x)
    ax2.set_xlabel("x")
    ax2.set_ylabel("ρ(x)")
    ax2.set_title("Densitas partikel tunggal (ground state)")
    st.pyplot(fig2)

    # Ringkasan
    st.markdown(f"""
    **Parameter simulasi:**
    - Panjang sumur (L) = {L}
    - Jumlah basis (Nmax) = {Nmax}
    - Kekuatan interaksi (g) = {g}
    - Grid resolusi (nx) = {nx}
    - Energi keadaan dasar = {Evals[0]:.6f}
    """)
else:
    st.info("Klik tombol **Jalankan Simulasi** untuk memulai perhitungan.")
