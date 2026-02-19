import streamlit as st
import requests
import re
from Bio import Entrez, SeqIO
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
fetch_button = st.button("Fetch & Analyze", on_click=trigger_analysis)

# -------------------------
# Run analysis
# -------------------------
if st.session_state.run_analysis:
    st.session_state.run_analysis = False

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
            st.success("Sequence retrieved successfully!")
            st.write(f"Organism: {organism}")
            st.write(f"Length: {len(sequence)} aa")

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
            # Extinction coefficient
            # -------------------------
            try:
                ec_reduced, ec_disulfide = analysis.molar_extinction_coefficient()
                st.subheader("Extinction Coefficient")
                st.write(f"Reduced cysteines: {ec_reduced:,} M⁻¹·cm⁻¹")
                st.write(f"Disulfide cysteines: {ec_disulfide:,} M⁻¹·cm⁻¹")
            except Exception as e:
                st.warning(f"Could not calculate extinction coefficient: {e}")

            # -------------------------
            # CYSTEINE-ANCHORED MOTIF DETECTION (all clusters)
            # -------------------------
            st.subheader("Detected Fe-S & Rubredoxin Motifs")

            cys_positions = [i for i, aa in enumerate(sequence) if aa == "C"]
            detected_regions = []

            def classify_cluster_by_spacings(cys_spacings):
                # Unpack spacings
                n = len(cys_spacings)
                if n == 3:
                    s1, s2, s3 = cys_spacings

                    # Classic Fe-S clusters
                    if 2 <= s1 <= 12 and 2 <= s2 <= 12 and 2 <= s3 <= 12:
                        return "4Fe-4S"
                    if 1 <= s1 <= 10 and 20 <= s2 <= 60 and 1 <= s3 <= 10:
                        return "2Fe-2S"
                    if 2 <= s1 <= 4 and 15 <= s2 <= 40 and 2 <= s3 <= 4:
                        return "Rubredoxin (long)"
                    if 1 <= s1 <= 4 and 4 <= s2 <= 12 and 1 <= s3 <= 4:
                        return "Rubredoxin (short)"
                    if s1 == 3 and s2 == 2 and s3 == 2:
                        return "Radical SAM (CxxxCxxC)"
                    if 1 <= s1 <= 3 and 15 <= s2 <= 30 and 1 <= s3 <= 3:
                        return "Rieske [2Fe-2S] (Cys-His)"
                    if 2 <= s1 <= 2 and 20 <= s2 <= 35 and s3 >= 1:
                        return "FeFe Hydrogenase H-cluster"
                elif n == 4:
                    s1, s2, s3, s4 = cys_spacings
                    # A-cluster
                    if 10 <= s1 <= 40 and 2 <= s2 <= 2 and 10 <= s3 <= 40:
                        return "A-cluster (ACS candidate)"
                    # C-cluster
                    if all(5 <= s <= 25 for s in (s1, s2, s3, s4)):
                        return "C-cluster (CODH candidate)"
                elif n == 5:
                    # P-cluster (Nitrogenase)
                    if all(10 <= s <= 60 for s in cys_spacings):
                        return "P-cluster (Nitrogenase candidate)"
                return None

            # Loop through cysteines as anchors
            for i in range(len(cys_positions)):
                for j in range(3, 6):  # Check next 3–5 cysteines
                    if i + j >= len(cys_positions):
                        continue
                    selected_cys = cys_positions[i:i+j+1]
                    spacings = [selected_cys[k+1]-selected_cys[k]-1 for k in range(len(selected_cys)-1)]
                    cluster_type = classify_cluster_by_spacings(spacings)
                    if cluster_type:
                        detected_regions.append({
                            "type": cluster_type,
                            "start": selected_cys[0],
                            "end": selected_cys[-1]+1,
                            "motif": sequence[selected_cys[0]:selected_cys[-1]+1],
                            "spacings": spacings
                        })

            # Remove duplicates
            unique_regions = []
            for region in detected_regions:
                if not any(r["start"] == region["start"] and r["end"] == region["end"]
                           for r in unique_regions):
                    unique_regions.append(region)

            if unique_regions:
                for idx, region in enumerate(unique_regions, 1):
                    st.markdown(f"### Cluster {idx}")
                    st.write(f"Type: {region['type']}")
                    st.write(f"Position: {region['start']+1}-{region['end']}")
                    st.write(f"Spacing: {region['spacings']}")
                    st.code(region["motif"])
            else:
                st.write("No Fe-S or rubredoxin motifs detected.")

            # -------------------------
            # Cysteine density plot
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
                yaxis_title="Cysteine presence",
                yaxis=dict(range=[0, 1.1])
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

            st.markdown(
                f'<div style="overflow-x:auto; white-space:pre-wrap; font-family:monospace">{html_seq}</div>',
                unsafe_allow_html=True
            )
