# SIMBA – Selection Inference and Mapping for Biological Analysis

**SIMBA** is a codon‑aware, end‑to‑end Python package for detecting, quantifying, and visualizing selection pressures in viral genomes.  
It bridges the gap between statistical selection detection and biological storytelling by integrating:

-  **Population‑level codon mapping** using dN/dS estimation  
-  **Annotation overlays** for functional and structural domains  
-  **Robust error handling** for ambiguous bases, frame shifts, and alignment quirks  
-  **Publication‑ready outputs**: tables, static plots, and annotated graphics  
-  **CLI‑driven workflow** for reproducible, rapid analysis  

---

## 🚀 Quick Start

### 1. Installation
```bash
git clone https://github.com/<YourUsername>/SIMBA.git
cd SIMBA
pip install -e .

python -m simba.cli \
  --fasta data/example.fasta \
  --annotation data/spike_annotations.csv \
  --output results/

SIMBA/
├── simba/                  # Core package modules
├── data/                   # Sample datasets & annotations
├── results/                # Output folder (user‑generated)
├── setup.py
├── README.md
└── .gitignore



---

## 📄 Create `.gitignore` in the same root folder
Save this as **`.gitignore`** right alongside your `README.md` and `setup.py`.

```gitignore
# Byte-compiled / cache
__pycache__/
*.py[cod]
*.so

# Virtual environments
.venv/
env/
venv/
ENV/

# OS-specific
.DS_Store
Thumbs.db
desktop.ini

# IDE / editor
.vscode/
.idea/

# Logs & temp
*.log
*.tmp

# Large raw data
*.fasta
*.fa
*.fastq
*.fq
*.bam
*.sam
*.vcf
*.gz
*.zip

# Generated results
results/
output/

# Jupyter checkpoints
.ipynb_checkpoints/
