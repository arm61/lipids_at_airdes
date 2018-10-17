current_dir = $(shell pwd)
FIG_DIR = reports/figures
PLOT_DIR = src/visualization
DLPC_DIR = output/dlpc
DMPC_DIR = output/dmpc
DPPC_DIR = output/dppc
DMPG_DIR = output/dmpg
DATA_DIR = data/processed

DLPC_DATA = $(DATA_DIR)/dlpc_xrr_sp_20.dat $(DATA_DIR)/dlpc_xrr_sp_25.dat $(DATA_DIR)/dlpc_xrr_sp_30.dat $(DATA_DIR)/dlpc_xrr_sp_35.dat
DMPC_DATA = $(DATA_DIR)/dmpc_xrr_sp_20.dat $(DATA_DIR)/dmpc_xrr_sp_25.dat $(DATA_DIR)/dmpc_xrr_sp_30.dat $(DATA_DIR)/dmpc_xrr_sp_40.dat
DPPC_DATA = $(DATA_DIR)/dppc_xrr_sp_15.dat $(DATA_DIR)/dppc_xrr_sp_20.dat $(DATA_DIR)/dppc_xrr_sp_25.dat $(DATA_DIR)/dlpc_xrr_sp_30.dat
DMPG_DATA = $(DATA_DIR)/dmpg_xrr_sp_15.dat $(DATA_DIR)/dmpg_xrr_sp_20.dat $(DATA_DIR)/dmpg_xrr_sp_25.dat $(DATA_DIR)/dlpc_xrr_sp_30.dat

NDMPC2_DATA = $(DATA_DIR)/dmpc_nr_hd_sp_20.dat $(DATA_DIR)/dmpc_nr_h_sp_20.dat
NDMPC1_DATA = $(DATA_DIR)/dmpc_nr_hd_sp_25.dat $(DATA_DIR)/dmpc_nr_h_sp_25.dat
NDPPC2_DATA = $(DATA_DIR)/dppc_nr_hd_sp_15.dat $(DATA_DIR)/dppc_nr_h_sp_15.dat
NDPPC1_DATA = $(DATA_DIR)/dppc_nr_hd_sp_20.dat $(DATA_DIR)/dppc_nr_h_sp_20.dat

PAPER_FIGS = $(FIG_DIR)/head_groups.png $(FIG_DIR)/dlpc_ref_sld.pdf $(FIG_DIR)/dmpc_ref_sld.pdf $(FIG_DIR)/dppc_ref_sld.pdf $(FIG_DIR)/dmpg_ref_sld.pdf $(FIG_DIR)/dlpc_vh_dt_phi.pdf $(FIG_DIR)/dmpc_vh_dt_phi.pdf $(FIG_DIR)/dppc_vh_dt_phi.pdf $(FIG_DIR)/dmpg_vh_dt_phi.pdf $(FIG_DIR)/dmpc_20n_ref_sld.pdf $(FIG_DIR)/dppc_20n_ref_sld.pdf

DLPC_PDF_FIGS = $(FIG_DIR)/dlpc1_all_corner.pdf $(FIG_DIR)/dlpc2_all_corner.pdf $(FIG_DIR)/dlpc3_all_corner.pdf $(FIG_DIR)/dlpc4_all_corner.pdf
DMPC_PDF_FIGS = $(FIG_DIR)/dmpc1_all_corner.pdf $(FIG_DIR)/dmpc2_all_corner.pdf $(FIG_DIR)/dmpc3_all_corner.pdf $(FIG_DIR)/dmpc4_all_corner.pdf
DPPC_PDF_FIGS = $(FIG_DIR)/dppc1_all_corner.pdf $(FIG_DIR)/dppc2_all_corner.pdf $(FIG_DIR)/dppc3_all_corner.pdf $(FIG_DIR)/dppc4_all_corner.pdf
DMPG_PDF_FIGS = $(FIG_DIR)/dmpg1_all_corner.pdf $(FIG_DIR)/dmpg2_all_corner.pdf $(FIG_DIR)/dmpg3_all_corner.pdf $(FIG_DIR)/dmpg4_all_corner.pdf
NDPPC_PDF_FIGS = $(FIG_DIR)/dppc_15n_all_corner.pdf $(FIG_DIR)/dppc_20n_all_corner.pdf
NDMPC_PDF_FIGS = $(FIG_DIR)/dmpc_20n_all_corner.pdf $(FIG_DIR)/dmpc_25n_all_corner.pdf
PDF_FIGS = $(DLPC_PDF_FIGS) $(DMPC_PDF_FIGS) $(DPPC_PDF_FIGS) $(DMPG_PDF_FIGS) $(NDPPC_PDF_FIGS) $(NDMPC_PDF_FIGS)

GIXD_FIGS = $(FIG_DIR)/gixd.png

ESI_FIGS = $(PDF_FIGS) $(GIXD_FIGS)


DLPC_OUT_PAPER = output/dlpc/angle35.txt output/dlpc/rough35.txt output/dlpc/vt.txt output/dlpc/vh.txt output/dlpc/head.txt output/dlpc/solh35.txt output/dlpc/tail35.txt
DMPC_OUT_PAPER = output/dmpc/angle40.txt output/dmpc/rough40.txt output/dmpc/vt.txt output/dmpc/vh.txt output/dmpc/head.txt output/dmpc/solh40.txt output/dmpc/tail40.txt output/dmpc/tail30.txt
DPPC_OUT_PAPER = output/dppc/angle30.txt output/dppc/rough30.txt output/dppc/vt.txt output/dppc/vh.txt output/dppc/head.txt output/dppc/solh30.txt output/dppc/tail30.txt
DMPG_OUT_PAPER = output/dmpg/angle30.txt output/dmpg/rough30.txt output/dmpg/vt.txt output/dmpg/vh.txt output/dmpg/head.txt output/dmpg/solh30.txt output/dmpg/tail30.txt

NDMPC1_OUT_PAPER = output/dmpc/angle25_neutron.txt output/dmpc/rough25_neutron.txt output/dmpc/solh25_neutron.txt output/dmpc/tail25_neutron.txt
NDMPC2_OUT_PAPER = output/dmpc/angle20_neutron.txt output/dmpc/rough20_neutron.txt output/dmpc/solh20_neutron.txt output/dmpc/tail20_neutron.txt
NDPPC1_OUT_PAPER = output/dppc/angle20_neutron.txt output/dppc/rough20_neutron.txt output/dppc/solh20_neutron.txt output/dppc/tail20_neutron.txt
NDPPC2_OUT_PAPER = output/dppc/angle15_neutron.txt output/dppc/rough15_neutron.txt output/dppc/solh15_neutron.txt output/dppc/tail15_neutron.txt

PAPER_OUT = $(DLPC_OUT_PAPER) $(DMPC_OUT_PAPER) $(DPPC_OUT_PAPER) $(DMPG_OUT_PAPER) $(NDMPC1_OUT_PAPER) $(NDMPC2_OUT_PAPER) $(NDPPC1_OUT_PAPER) $(NDPPC2_OUT_PAPER)


DLPC_OUT_ESI =output/dlpc/angle20.txt output/dlpc/rough20.txt output/dlpc/solh20.txt output/dlpc/tail20.txt output/dlpc/angle25.txt output/dlpc/rough25.txt  output/dlpc/solh25.txt output/dlpc/tail25.txt output/dlpc/angle30.txt output/dlpc/rough30.txt output/dlpc/solh30.txt output/dlpc/tail30.txt
DMPC_OUT_ESI =output/dmpc/angle20.txt output/dmpc/rough20.txt output/dmpc/solh20.txt output/dmpc/tail20.txt output/dmpc/angle25.txt output/dmpc/rough25.txt  output/dmpc/solh25.txt output/dmpc/tail25.txt output/dmpc/angle30.txt output/dmpc/rough30.txt output/dmpc/solh30.txt
DPPC_OUT_ESI = output/dppc/angle15.txt output/dppc/rough15.txt output/dppc/solh15.txt output/dppc/tail15.txt output/dppc/angle20.txt output/dppc/rough20.txt  output/dppc/solh20.txt output/dppc/tail20.txt output/dppc/angle25.txt output/dppc/rough25.txt output/dppc/solh25.txt output/dppc/tail25.txt
DMPG_OUT_ESI = output/dmpg/angle15.txt output/dmpg/rough15.txt output/dmpg/solh15.txt output/dmpg/tail15.txt output/dmpg/angle20.txt output/dmpg/rough20.txt  output/dmpg/solh20.txt output/dmpg/tail20.txt output/dmpg/angle25.txt output/dmpg/rough25.txt output/dmpg/solh25.txt output/dmpg/tail25.txt

ESI_OUT = $(DLPC_OUT_ESI) $(DMPC_OUT_ESI) $(DPPC_OUT_ESI) $(DMPG_OUT_ESI)


all : reports/si.pdf reports/paper.pdf
clean :
	rm -r reports/paper.pdf reports/si.pdf $(FIG_DIR)/*ref_sld.pdf  $(FIG_DIR)/*vh_dt_phi.pdf $(FIG_DIR)/*all_corner.pdf $(PLOT_DIR)/*.py notebooks/*.py output/dlpc/* output/dmpc/* output/dppc/* output/dmpg/* 

reports/si.pdf : reports/si.tex reports/bibi.bib $(ESI_FIGS) $(ESI_OUT)
	cd reports && pdflatex si.tex
	cd reports && bibtex si.aux
	cd reports && pdflatex si.tex
	cd reports && pdflatex si.tex
reports/paper.pdf : reports/paper.tex reports/bibi.bib $(PAPER_FIG) $(PAPER_OUT) reports/si.pdf
	cd reports && pdflatex paper.tex
	cd reports && bibtex paper.aux
	cd reports && pdflatex paper.tex
	cd reports && pdflatex paper.tex

output/dlpc/chain.txt : notebooks/lipid_xrr.py src/models/mol_vol.py $(DLPC_DATA)
	cd notebooks && ipython lipid_xrr.py dlpc 11 20 25 30 35 a
$(DLPC_OUT_PAPER) $(DLPC_OUT_ESI) $(FIG_DIR)/DLPC_ref_sld.pdf $(FIG_DIR)/DLPC_vh_dt_phi.pdf $(DLPC_PDF_FIGS): notebooks/xrr_chain_analysis.py src/models/mol_vol.py output/dlpc/chain.txt
	cd notebooks && ipython xrr_chain_analysis.py dlpc 11 20 25 30 35 a
output/dmpc/chain.txt : notebooks/lipid_xrr.py src/models/mol_vol.py $(DMPC_DATA)
	cd notebooks && ipython lipid_xrr.py dmpc 13 20 25 30 40 b
$(DMPC_OUT_PAPER) $(DMPC_OUT_ESI) $(FIG_DIR)/DMPC_ref_sld.pdf $(FIG_DIR)/DMPC_vh_dt_phi.pdf $(DMPC_PDF_FIGS): notebooks/xrr_chain_analysis.py src/models/mol_vol.py output/dmpc/chain.txt
	cd notebooks && ipython xrr_chain_analysis.py dmpc 13 20 25 30 40 b
output/dppc/chain.txt : notebooks/lipid_xrr.py src/models/mol_vol.py $(DPPC_DATA)
	cd notebooks && ipython lipid_xrr.py dppc 15 15 20 25 30 c
$(DPPC_OUT_PAPER) $(DPPC_OUT_ESI) $(FIG_DIR)/DPPC_ref_sld.pdf $(FIG_DIR)/DPPC_vh_dt_phi.pdf $(DPPC_PDF_FIGS): notebooks/xrr_chain_analysis.py src/models/mol_vol.py output/dppc/chain.txt
	cd notebooks && ipython xrr_chain_analysis.py dppc 15 15 20 25 30 c
output/dmpg/chain.txt : notebooks/lipid_xrr.py src/models/mol_vol.py $(DMPG_DATA)
	cd notebooks && ipython lipid_xrr.py dmpg 13 15 20 25 30 d
$(DMPG_OUT_PAPER) $(DMPG_OUT_ESI) $(FIG_DIR)/DMPG_ref_sld.pdf $(FIG_DIR)/DMPG_vh_dt_phi.pdf $(DMPG_PDF_FIGS): notebooks/xrr_chain_analysis.py src/models/mol_vol.py output/dmpg/chain.txt
	cd notebooks && ipython xrr_chain_analysis.py dmpg 13 15 20 25 30 d


output/dmpc/25_chain_neutron.txt : notebooks/lipid_nr.py src/models/mol_vol.py $(DMPC1_DATA)
	cd notebooks && ipython lipid_nr.py dmpc 13 25 a
$(NDMPC1_OUT_PAPER) $(FIG_DIR)/dmpc_25n_ref_sld.pdf $(FIG_DIR)/dmpc_25n_all_corner.pdf : notebooks/nr_chain_analysis.py src/models/mol_vol.py output/dmpc/25_chain_neutron.txt
	cd notebooks && ipython nr_chain_analysis.py dmpc 13 25 a
output/dmpc/20_chain_neutron.txt : notebooks/lipid_nr.py src/models/mol_vol.py $(DMPC2_DATA)
	cd notebooks && ipython lipid_nr.py dmpc 13 20 a
$(NDMPC2_OUT_PAPER) $(FIG_DIR)/dmpc_20n_ref_sld.pdf $(FIG_DIR)/dmpc_20n_all_corner.pdf : notebooks/nr_chain_analysis.py src/models/mol_vol.py output/dmpc/20_chain_neutron.txt
	cd notebooks && ipython nr_chain_analysis.py dmpc 13 20 a
output/dppc/20_chain_neutron.txt : notebooks/lipid_nr.py src/models/mol_vol.py $(DPPC1_DATA)
	cd notebooks && ipython lipid_nr.py dppc 15 20 b
$(NDPPC1_OUT_PAPER) $(FIG_DIR)/dppc_20n_ref_sld.pdf $(FIG_DIR)/dppc_20n_all_corner.pdf : notebooks/nr_chain_analysis.py src/models/mol_vol.py output/dppc/20_chain_neutron.txt
	cd notebooks && ipython nr_chain_analysis.py dppc 15 20 b
output/dppc/15_chain_neutron.txt : notebooks/lipid_nr.py src/models/mol_vol.py $(DPPC2_DATA)
	cd notebooks && ipython lipid_nr.py dppc 15 15 b
$(NDPPC2_OUT_PAPER) $(FIG_DIR)/dppc_15n_ref_sld.pdf $(FIG_DIR)/dppc_15n_all_corner.pdf : notebooks/nr_chain_analysis.py src/models/mol_vol.py output/dppc/15_chain_neutron.txt
	cd notebooks && ipython nr_chain_analysis.py dppc 15 15 b


notebooks/lipid_xrr.py : notebooks/lipid_xrr.ipynb
	jupyter-nbconvert notebooks/lipid_xrr.ipynb --to script --output-dir=notebooks/
notebooks/lipid_nr.py : notebooks/lipid_nr.ipynb
	jupyter-nbconvert notebooks/lipid_nr.ipynb --to script --output-dir=notebooks/
notebooks/xrr_chain_analysis.py : notebooks/xrr_chain_analysis.ipynb
	jupyter-nbconvert notebooks/xrr_chain_analysis.ipynb --to script --output-dir=notebooks/
notebooks/nr_chain_analysis.py : notebooks/nr_chain_analysis.ipynb
	jupyter-nbconvert notebooks/nr_chain_analysis.ipynb --to script --output-dir=notebooks/
