import streamlit as st
from Bio import Entrez, SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

Entrez.email = "nima.ninate@gmail.com"

st.title("🧬 GeneCompare: Cross-Species Phylogenetics Explorer")

# Initialize session state
if "show_results" not in st.session_state:
    st.session_state["show_results"] = False

st.sidebar.header("Gene Input")
gene = st.sidebar.text_input("Gene name (e.g., BRCA1, COI)", "COI")
species_input = st.sidebar.text_area("List up to 10 species (one per line)",
                                     "Homo sapiens\nPan troglodytes\nFelis catus")
species_list = [s.strip() for s in species_input.strip().splitlines() if s.strip()][:10]

if st.sidebar.button("Fetch & Analyze"):
    st.header(f"Analysis for gene: {gene}")

    os.makedirs("data", exist_ok=True)

    sequences = []
    errors = []
    for sp in species_list:
        term = f"{gene}[Gene] AND {sp}[Organism]"
        try:
            handle = Entrez.esearch(db="nucleotide", term=term, retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                errors.append(f"No results for {sp}")
                continue

            id_ = record["IdList"][0]
            handle = Entrez.efetch(db="nucleotide", id=id_, rettype="fasta", retmode="text")
            seq_record = SeqIO.read(handle, "fasta")
            handle.close()

            seq_record.id = sp.replace(" ", "_")
            sequences.append(seq_record)
        except Exception as e:
            errors.append(f"Error fetching {sp}: {str(e)}")

    if errors:
        st.error("\n".join(errors))

    if len(sequences) < 2:
        st.warning("Not enough sequences to align.")
        st.session_state["show_results"] = False
    else:
        # Translate DNA to Protein
        for i, rec in enumerate(sequences):
            dna_seq = str(rec.seq).replace("\n", "").replace(" ", "")
            dna_seq = dna_seq.upper().replace("T", "U")  # RNA style
            start = dna_seq.find("AUG")
            protein = ""
            codon_table = {
                codon: str(Seq(codon).transcribe().translate())
                for codon in [a + b + c for a in "UCAG" for b in "UCAG" for c in "UCAG"]
            }
            for j in range(start, len(dna_seq), 3):
                codon = dna_seq[j:j + 3]
                if len(codon) != 3:
                    break
                aa = codon_table.get(codon, "")
                if aa == "*":
                    break
                protein += aa
            rec.seq = Seq(protein)

        # Temporary workaround → force all sequences to same length
        min_len = min(len(rec.seq) for rec in sequences)
        trimmed_sequences = [
            SeqRecord(rec.seq[:min_len], id=rec.id, description="")
            for rec in sequences
        ]
        alignment = MultipleSeqAlignment(trimmed_sequences)

        # Save alignment to file
        with open("data/alignment.fasta", "w") as f:
            SeqIO.write(alignment, f, "fasta")

        # Distance matrix
        calculator = DistanceCalculator("blosum62")
        dm = calculator.get_distance(alignment)
        df = pd.DataFrame(dm.matrix, index=dm.names, columns=dm.names)

        # Save data to session state for stable downloads and stable plots
        st.session_state["alignment_data"] = open("data/alignment.fasta").read()
        st.session_state["identity_csv"] = df.to_csv()
        st.session_state["distance_matrix"] = df
        st.session_state["tree"] = DistanceTreeConstructor().nj(dm)
        st.session_state["show_results"] = True

# Show results if available
if st.session_state["show_results"]:
    st.subheader("Pairwise Identity Heatmap")
    fig, ax = plt.subplots()
    sns.heatmap(1 - st.session_state["distance_matrix"], annot=True, cmap="viridis", ax=ax)
    st.pyplot(fig)

    st.subheader("Phylogenetic Tree")
    fig2, ax2 = plt.subplots(figsize=(8, 6))
    Phylo.draw(st.session_state["tree"], do_show=False, axes=ax2)
    st.pyplot(fig2)

# Download buttons (always visible if data present)
if "alignment_data" in st.session_state:
    st.download_button("Download Alignment (FASTA)", data=st.session_state["alignment_data"], file_name="alignment.fasta")

if "identity_csv" in st.session_state:
    st.download_button("Download Identity Matrix (CSV)", data=st.session_state["identity_csv"], file_name="identity_matrix.csv")
