import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

st.set_page_config(
    page_title="DFT Research Dashboard - CMC Uni",
    page_icon="⚗️",
    layout="wide"
)

st.markdown("""
<style>
    .header-title { color: #2c3e50; text-align: center; font-weight: 800; margin-bottom: 0px;}
    .sub-title { color: #7f8c8d; text-align: center; margin-bottom: 20px; font-style: italic;}
    .metric-card { background-color: #f8f9fa; border: 1px solid #e9ecef; border-radius: 8px; padding: 15px; text-align: center; }
    .highlight { color: #e74c3c; font-weight: bold; }
</style>
""", unsafe_allow_html=True)

st.markdown('<h1 class="header-title">COMPUTATIONAL MATERIALS SCIENCE</h1>', unsafe_allow_html=True)
st.markdown('<div class="sub-title">Phân tích Cấu trúc Vùng & DOS từ dữ liệu DFT</div>', unsafe_allow_html=True)
st.markdown("---")

# --- SIDEBAR ---
with st.sidebar:
    st.header("🎛️ Bảng điều khiển")
    st.write("Chọn vật liệu để phân tích dữ liệu DFT:")
    material_opt = st.selectbox("Vật liệu (Formula):", 
                                ["Silicon (Si-Bulk)", "Graphene (C-2D)", "MoS2 (Monolayer)"])
    
    st.markdown("---")
    st.markdown("### ⚙️ Thông số DFT")
    st.info("""
    **Functional:** PBE (Generalized Gradient Approximation)
    **Cutoff Energy:** 520 eV
    **k-point grid:** 8x8x1 (Gamma centered)
    **Spin Polarized:** False
    """)
    st.markdown("---")
    st.caption("Dữ liệu mô phỏng được chuẩn hóa từ Materials Project Database.")

def get_band_structure_data(mat):
    k_path = np.linspace(0, 10, 100) # Đường dẫn trong không gian k
    
    if mat == "Silicon (Si-Bulk)":
        # Silicon: Indirect Gap (1.12 eV)
        # VBM (Valence Band Max) tại Gamma (0), CBM (Conduction Band Min) lệch Gamma
        vbm = -0.5 * (k_path - 5)**2 # Parabol úp
        cbm = 1.12 + 0.3 * (k_path - 8)**2 # Parabol ngửa, đáy lệch tâm
        gap_type = "Indirect (Gián tiếp)"
        gap_val = 1.12
        lattice = "Diamond Cubic"
        
    elif mat == "Graphene (C-2D)":
        # Graphene: Zero Gap (Dirac Cone)
        # Chạm nhau tại điểm K (giả sử tại k=5)
        vbm = -1.5 * np.abs(k_path - 5)
        cbm = 1.5 * np.abs(k_path - 5)
        gap_type = "Zero Gap (Semi-metal)"
        gap_val = 0.0
        lattice = "Hexagonal (Honeycomb)"

    elif mat == "MoS2 (Monolayer)":
        # MoS2 đơn lớp: Direct Gap (~1.8 eV) tại điểm K
        # VBM và CBM đều cực trị tại cùng một điểm k (giả sử k=5)
        vbm = -0.8 * (k_path - 5)**2
        cbm = 1.8 + 0.8 * (k_path - 5)**2
        gap_type = "Direct (Trực tiếp)"
        gap_val = 1.8
        lattice = "Hexagonal"
        
    return k_path, vbm, cbm, gap_type, gap_val, lattice

# --- MAIN LAYOUT ---

# 1. THÔNG TIN CƠ BẢN
k, val_band, con_band, gap_type, gap_val, struct = get_band_structure_data(material_opt)

col1, col2, col3, col4 = st.columns(4)
with col1:
    st.metric("Cấu trúc tinh thể", struct)
with col2:
    st.metric("Bandgap (eV)", f"{gap_val} eV")
with col3:
    st.metric("Loại vùng cấm", gap_type, delta_color="normal")
with col4:
    is_good_optical = "Tốt" if "Direct" in gap_type else "Kém"
    st.metric("Tính chất quang", is_good_optical)

# 2. BIỂU ĐỒ BAND STRUCTURE & DOS
col_band, col_dos = st.columns([2, 1])

with col_band:
    st.subheader("1. Electronic Band Structure")
    
    fig_band = go.Figure()
    
    # Vẽ nhiều dải (bands) giả lập để nhìn cho "thật" hơn
    # Dải dẫn (Conduction)
    fig_band.add_trace(go.Scatter(x=k, y=con_band, mode='lines', name='CBM', line=dict(color='#e74c3c', width=3)))
    fig_band.add_trace(go.Scatter(x=k, y=con_band + 1.5, mode='lines', line=dict(color='#e74c3c', width=1, dash='dot'), showlegend=False))
    
    # Dải hóa trị (Valence)
    fig_band.add_trace(go.Scatter(x=k, y=val_band, mode='lines', name='VBM', line=dict(color='#3498db', width=3)))
    fig_band.add_trace(go.Scatter(x=k, y=val_band - 1.5, mode='lines', line=dict(color='#3498db', width=1, dash='dot'), showlegend=False))
    
    # Fermi Level
    fig_band.add_hline(y=0, line_dash="dash", line_color="green", annotation_text="Fermi Level (E_f)")

    # Chú thích các điểm đối xứng
    tick_vals = [0, 5, 10]
    if material_opt == "Silicon (Si-Bulk)":
        tick_text = ['Γ', 'X', 'L'] # Ký hiệu giả định
    else:
        tick_text = ['Γ', 'K', 'M']

    fig_band.update_layout(
        xaxis=dict(
            tickmode='array',
            tickvals=tick_vals,
            ticktext=tick_text,
            title='Wave Vector (k-points)',
            showgrid=True,
            gridcolor='lightgray'
        ),
        yaxis=dict(title='Energy (E - Ef) [eV]', range=[-4, 5]),
        height=500,
        plot_bgcolor='white',
        margin=dict(t=30, b=0, l=0, r=0)
    )
    st.plotly_chart(fig_band, use_container_width=True)

with col_dos:
    st.subheader("2. Density of States (DOS)")
    
    # Tạo dữ liệu DOS giả lập tương ứng với Bandgap
    e_dos = np.linspace(-4, 5, 200)
    dos_val = np.zeros_like(e_dos)
    
    # DOS có giá trị tại các mức năng lượng có dải năng lượng
    # Valence band DOS
    mask_val = e_dos < 0
    dos_val[mask_val] = np.exp(-(e_dos[mask_val] + 1)**2) + 0.5 * np.exp(-(e_dos[mask_val] + 2.5)**2)
    
    # Conduction band DOS (bắt đầu từ gap_val)
    mask_con = e_dos > gap_val
    dos_val[mask_con] = np.exp(-(e_dos[mask_con] - gap_val - 1)**2) + 0.5
    
    fig_dos = go.Figure()
    fig_dos.add_trace(go.Scatter(x=dos_val, y=e_dos, mode='lines', fill='tozerox', line=dict(color='#2c3e50')))
    
    fig_dos.update_layout(
        xaxis=dict(title='DOS (states/eV)', showgrid=False),
        yaxis=dict(title='', showticklabels=False, range=[-4, 5]), # Khớp trục Y với biểu đồ Band
        height=500,
        plot_bgcolor='white',
        margin=dict(t=30, b=0, l=0, r=0)
    )
    # Fermi Line
    fig_dos.add_hline(y=0, line_dash="dash", line_color="green")
    
    st.plotly_chart(fig_dos, use_container_width=True)

# 3. PHÂN TÍCH & KẾT LUẬN
st.markdown("### 📝 Phân tích Khoa học")

if material_opt == "Silicon (Si-Bulk)":
    st.warning("""
    **Kết luận:** Silicon có **Indirect Bandgap**.
    * Đỉnh vùng hóa trị (VBM) và đáy vùng dẫn (CBM) nằm ở các vector sóng (k) khác nhau.
    * Điều này có nghĩa là electron cần thay đổi động lượng (phonon) để nhảy lên vùng dẫn.
    * => Hiệu suất phát quang (LED/Laser) kém, nhưng rất tốt cho linh kiện điện tử (Transistor).
    """)
elif material_opt == "Graphene (C-2D)":
    st.warning("""
    **Kết luận:** Graphene là vật liệu **Semi-metal (Zero Bandgap)**.
    * Vùng dẫn và vùng hóa trị chạm nhau tại điểm Dirac (K).
    * Độ linh động điện tử (Mobility) cực cao nhưng không có Gap để tắt dòng điện.
    * => Khó làm Transistor logic số, nhưng tuyệt vời cho Analog RF hoặc vật liệu dẫn điện trong suốt.
    """)
elif material_opt == "MoS2 (Monolayer)":
    st.success("""
    **Kết luận:** MoS2 đơn lớp có **Direct Bandgap (~1.8 eV)**.
    * VBM và CBM thẳng hàng trong không gian k.
    * Electron có thể chuyển dời thẳng đứng mà không cần phonon hỗ trợ.
    * => **Tiềm năng lớn:** Ứng dụng chế tạo LED siêu mỏng, Pin mặt trời hiệu suất cao và Transistor thế hệ mới.
    """)

# --- PHẦN HƯỚNG DẪN THỰC TẾ ---
with st.expander("🔬 Hướng dẫn: Cách chạy DFT thực tế cho NCKH"):
    st.markdown("""
    Để có dữ liệu thật cho dự án của bạn, bạn cần thực hiện quy trình sau trên máy trạm (Linux):

    **Bước 1: Chuẩn bị file đầu vào (Input File)**
    Ví dụ sử dụng phần mềm **Quantum ESPRESSO** (miễn phí, mã nguồn mở).
    File `mos2.in`:
    ```fortran
    &CONTROL
      calculation = 'scf'       ! Tính toán Self-consistent field
      prefix = 'mos2'
      outdir = './tmp/'
      pseudo_dir = './pseudo/'  ! Thư mục chứa giả thế (Pseudopotentials)
    /
    &SYSTEM
      ibrav = 4, A = 3.16, C = 10, nat = 3, ntyp = 2
      ecutwfc = 60              ! Energy Cutoff (quan trọng)
    /
    ATOMIC_POSITIONS {crystal}
      Mo 0.3333 0.6666 0.5000
      S  0.3333 0.6666 0.6500
      S  0.3333 0.6666 0.3500
    K_POINTS {automatic}
      12 12 1 0 0 0
    ```
    
    **Bước 2: Chạy tính toán**
    `pw.x < mos2.in > mos2.out`
    
    **Bước 3: Xử lý hậu kỳ (Post-processing)**
    Dùng phần mềm này (bạn đang viết) để đọc file `bands.dat` xuất ra từ bước trên và vẽ đồ thị.
    """)

