import streamlit as st
import requests
import re
from Bio import Entrez, SeqIO, __version__ as biopython_version
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import plotly.graph_objects as go

st.set_page_config(page_title="Fe-S & Rubredoxin Analyzer", layout="wide")
st.title("Fe-S & Rubredoxin Analyzer")

# -------------------------
# Session state for Enter key
# -------------------------
if "run_analysis" not in st.session_state:
    st.session_state.run_analysis = False

def trigger_analysis():
    st.session_state.run_analysis = True

# -------------------------
# Input selection
# -------------------------
source = st.selectbox("Select input source", ["UniProt", "NCBI", "KEGG"])
protein_id = st.text_input("Enter protein ID", key="protein_input", on_change=trigger_analysis)

# Fetch & Analyze button
fetch_button = st.button("Fetch & Analyze", on_click=trigger_analysis)

# Run analysis if triggered
if st.session_state.run_analysis:
    st.session_state.run_analysis = False  # reset trigger

    if not protein_id:
        st.warning("Please enter a protein ID.")
    else:
        sequence = ""
        header = ""
        organism = "Unknown"

        # -------------------------
        # Fetch sequence
        # -------------------------
        try:
            if source == "UniProt":
                url = f"https://rest.uniprot.org/uniprotkb/{protein_id}.fasta"
                response = requests.get(url)
                response.raise_for_status()
                fasta = response.text
                header = fasta.split("\n")[0]
                sequence = "".join(fasta.split("\n")[1:]).upper()
                match = re.search(r'OS=(.*?)\s(?:GN=|PE=|SV=|$)', header)
                if match:
                    organism = match.group(1)

            elif source == "NCBI":
                Entrez.email = "your_email@example.com"
                handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
                record = SeqIO.read(handle, "fasta")
                handle.close()
                sequence = str(record.seq).upper()
                header = record.description
                match = re.search(r'\[(.*?)\]', header)
                if match:
                    organism = match.group(1)

            elif source == "KEGG":
                url = f"http://rest.kegg.jp/get/{protein_id}/aaseq"
                response = requests.get(url)
                response.raise_for_status()
                lines = response.text.split("\n")
                header = lines[0]
                sequence = "".join(lines[1:]).upper()
                if ":" in protein_id:
                    organism = protein_id.split(":")[0]

        except Exception as e:
            st.error(f"Error fetching sequence: {e}")

        if sequence:
            st.success(f"Sequence retrieved successfully! Length: {len(sequence)} aa")
            st.write(f"Organism: {organism}")
            st.write(f"Biopython version: {biopython_version}")

            # -------------------------
            # Basic properties
            # -------------------------
            analysis = ProteinAnalysis(sequence)
            mw = analysis.molecular_weight()
            pI = analysis.isoelectric_point()
            cysteine_count = sequence.count("C")

            st.subheader("Basic Properties")
            st.write(f"Molecular Weight: {mw:,.2f} Da")
            st.write(f"Isoelectric Point (pI): {pI:.2f}")
            st.write(f"Total Cysteines: {cysteine_count}")

            # -------------------------
            # Extinction coefficient (tuple-based for Biopython 1.85)
            # -------------------------
            try:
                ec_reduced, ec_disulfide = analysis.molar_extinction_coefficient()
                st.subheader("Extinction Coefficient")
                st.write(f"Reduced cysteines: {ec_reduced:,} M⁻¹·cm⁻¹")
                st.write(f"Disulfide cysteines: {ec_disulfide:,} M⁻¹·cm⁻¹")
            except Exception as e:
                st.warning(f"Could not calculate extinction coefficient: {e}")

            # -------------------------
            # Hierarchical motif detection
            # -------------------------
            motif_definitions = [
                {"type": "4Fe-4S", "pattern": r"C.{2,12}C.{2,12}C.{2,20}C", "c_count": 4},
                {"type": "3Fe-4S", "pattern": r"C.{2,15}C.{2,15}C", "c_count": 3},
                {"type": "2Fe-2S", "pattern": r"C.{1,10}C.{1,10}C", "c_count": 2},
                {"type": "rubredoxin_short", "pattern": r"C.{1,4}C.{4,12}C.{1,4}C", "c_count": 4},
                {"type": "rubredoxin_long", "pattern": r"C.{2}C.{15,40}C.{2}C", "c_count": 4},
            ]

            assigned_motifs = []
            occupied_positions = set()

            for mf_def in motif_definitions:
                pattern = re.compile(mf_def["pattern"])
                for match in pattern.finditer(sequence):
                    start = match.start()
                    end = match.end()
                    if any(pos in occupied_positions for pos in range(start, end)):
                        continue
                    motif_seq = match.group()
                    c_positions = [m.start() for m in re.finditer("C", motif_seq)]
                    if len(c_positions) != mf_def["c_count"]:
                        continue
                    spacings = [c_positions[i+1]-c_positions[i]-1 for i in range(len(c_positions)-1)]
                    assigned_motifs.append({
                        "type": mf_def["type"],
                        "start": start,
                        "motif": motif_seq,
                        "spacings": spacings
                    })
                    occupied_positions.update(range(start, end))

            # -------------------------
            # Display motifs
            # -------------------------
            st.subheader("Detected Fe-S & Rubredoxin Motifs")
            if assigned_motifs:
                for mf in assigned_motifs:
                    st.write(f"{mf['type']} | Position {mf['start']+1}: {mf['motif']} | Spacings: {mf['spacings']}")
            else:
                st.write("No Fe-S or rubredoxin motifs detected.")

            # -------------------------
            # Cysteine density
            # -------------------------
            positions = list(range(1, len(sequence)+1))
            densities = [1 if aa=="C" else 0 for aa in sequence]

            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=positions,
                y=densities,
                mode="lines+markers",
                name="Cys presence",
                line=dict(color='black'),
                marker=dict(size=6)
            ))

            fig.update_layout(
                height=250,
                margin=dict(l=20, r=20, t=20, b=20),
                xaxis_title="Residue position",
                yaxis_title="Cysteine density",
                yaxis=dict(range=[0, max(densities)+1])
            )

            st.plotly_chart(fig)

            # -------------------------
            # Sequence highlighting
            # -------------------------
            st.subheader("Sequence (Cysteines highlighted in yellow)")
            html_seq = ""
            for aa in sequence:
                if aa == "C":
                    html_seq += f'<span style="background-color:yellow">{aa}</span>'
                else:
                    html_seq += aa
            st.markdown(f'<div style="overflow-x:auto; white-space:pre">{html_seq}</div>', unsafe_allow_html=True)
