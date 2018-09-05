current_dir = $(shell pwd)
FIG_DIR = reports/figures
PLOT_DIR = src/visualization
DLPC_DIR = output/dlpc
DMPC_DIR = output/dmpc
DPPC_DIR = output/dppc
DMPG_DIR = output/dmpg
DATA_DIR = data/processed

PAPER_FIGS = $(FIG_DIR)/dlpc_ref_sld.pdf $(FIG_DIR)/dmpc_ref_sld.pdf $(FIG_DIR)/dmpg_ref_sld.pdf $(FIG_DIR)/dppc_ref_sld.pdf  $(FIG_DIR)/dlpc_vh_dt_phi.pdf $(FIG_DIR)/dmpc_vh_dt_phi.pdf $(FIG_DIR)/dmpg_vh_dt_phi.pdf $(FIG_DIR)/dppc_vh_dt_phi.pdf $(FIG_DIR)/nDPPC15_all_data.pdf $(FIG_DIR)/nDMPC20_all_data.pdf $(FIG_DIR)/nDPPC20_all_data.pdf $(FIG_DIR)/nDMPC25_all_data.pdf $(FIG_DIR)/head_groups.png
ESI_FIGS = $(FIG_DIR)/dlpc4_all_corner.pdf $(FIG_DIR)/dlpc5_all_corner.pdf $(FIG_DIR)/dmpc25_neutron_corner.pdf $(FIG_DIR)/dmpc20_neutron_corner.pdf $(FIG_DIR)/dmpc4_all_corner.pdf $(FIG_DIR)/dmpc5_all_corner.pdf $(FIG_DIR)/dmpg4_all_corner.pdf $(FIG_DIR)/dmpg5_all_corner.pdf $(FIG_DIR)/dppc20_neutron_corner.pdf $(FIG_DIR)/dppc15_neutron_corner.pdf $(FIG_DIR)/dppc4_all_corner.pdf $(FIG_DIR)/dppc5_all_corner.pdf
DLPC_OUT = output/dlpc/chain.txt
DMPC_OUT = output/dmpc/chain.txt
DPPC_OUT = output/dppc/chain.txt
DMPG_OUT = output/dmpg/chain.txt
NDMPC1_OUT = output/dmpc/25_chain_neutron.txt
NDMPC2_OUT = output/dmpc/20_chain_neutron.txt
NDPPC1_OUT = output/dppc/20_chain_neutron.txt
NDPPC2_OUT = output/dppc/15_chain_neutron.txt

all : reports/paper.pdf reports/esi.pdf reports/apssamp.pdf reports/esi2.pdf
clean :
	rm reports/paper.pdf reports/esi.pdf $(FIG_DIR)/*ref_sld.pdf  $(FIG_DIR)/*vh_dt_phi.pdf $(FIG_DIR)/*all_corner.pdf $(FIG_DIR)/*neutron_corner*.pdf $(PLOT_DIR)/*.py notebooks/DLPC/*.py notebooks/DMPC/*.py notebooks/DPPC/*.py notebooks/DMPG/*.py output/* output/dlpc/* output/dppc/* output/dmpc/* output/dmpg/*

reports/paper.pdf : reports/paper.tex reports/rsc.bib $(PAPER_FIGS) $(DLPC_OUT) $(DMPC_OUT) $(DPPC_OUT) $(DMPG_OUT)
	cd reports && pdflatex paper.tex
	cd reports && bibtex paper.aux
	cd reports && pdflatex paper.tex
	cd reports && pdflatex paper.tex
reports/esi.pdf : reports/esi.tex $(ESI_FIGS)
	cd reports && pdflatex esi.tex
	cd reports && bibtex esi.aux
	cd reports && pdflatex esi.tex
	cd reports && pdflatex esi.tex
reports/apssamp.pdf : reports/apssamp.tex reports/rsc.bib $(PAPER_FIGS) $(DLPC_OUT) $(DMPC_OUT) $(DPPC_OUT) $(DMPG_OUT)
	cd reports && pdflatex apssamp.tex
	cd reports && bibtex apssamp.aux
	cd reports && pdflatex apssamp.tex
	cd reports && pdflatex apssamp.tex
reports/esi2.pdf : reports/esi2.tex $(ESI_FIGS)
	cd reports && pdflatex esi2.tex
	cd reports && bibtex esi2.aux
	cd reports && pdflatex esi2.tex
	cd reports && pdflatex esi2.tex


$(FIG_DIR)/nDMPC25_all_data.pdf : $(PLOT_DIR)/plotrefn.py $(NDMPC1_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dmpc25_1_ref_neutron.txt ../../output/dmpc25_1_sld_neutron.txt ../../output/dmpc25_2_ref_neutron.txt ../../output/dmpc25_2_sld_neutron.txt ../../output/dmpc25_chain_neutron.txt a nDMPC25
$(FIG_DIR)/nDMPC20_all_data.pdf : $(PLOT_DIR)/plotrefn.py $(NDMPC2_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dmpc20_1_ref_neutron.txt ../../output/dmpc20_1_sld_neutron.txt ../../output/dmpc20_2_ref_neutron.txt ../../output/dmpc20_2_sld_neutron.txt ../../output/dmpc20_chain_neutron.txt a nDMPC20
$(FIG_DIR)/nDPPC20_all_data.pdf : $(PLOT_DIR)/plotrefn.py $(NDPPC1_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dppc20_1_ref_neutron.txt ../../output/dppc20_1_sld_neutron.txt ../../output/dppc20_2_ref_neutron.txt ../../output/dppc20_2_sld_neutron.txt ../../output/dppc20_chain_neutron.txt b nDPPC20
$(FIG_DIR)/nDPPC15_all_data.pdf : $(PLOT_DIR)/plotrefn.py $(NDPPC2_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dppc15_1_ref_neutron.txt ../../output/dppc15_1_sld_neutron.txt ../../output/dppc15_2_ref_neutron.txt ../../output/dppc15_2_sld_neutron.txt ../../output/dppc15_chain_neutron.txt b nDPPC15

$(PLOT_DIR)/plotrefn.py : $(PLOT_DIR)/plotrefn.ipynb
	jupyter-nbconvert src/visualization/plotrefn.ipynb --to script --output-dir=src/visualization/

$(DLPC_OUT) $(FIG_DIR)/DLPC_ref_sld.pdf $(FIG_DIR)/DLPC_vh_dt_phi.pdf : notebooks/lipid_xrr.py src/models/mol_vol.py
	cd notebooks && ipython lipid_xrr.py dlpc 11 20 25 30 35 a
$(DMPC_OUT) $(FIG_DIR)/DMPC_ref_sld.pdf $(FIG_DIR)/DMPC_vh_dt_phi.pdf : notebooks/lipid_xrr.py src/models/mol_vol.py
	cd notebooks && ipython lipid_xrr.py dmpc 13 20 25 30 40 b
$(DPPC_OUT) $(FIG_DIR)/DPPC_ref_sld.pdf $(FIG_DIR)/DPPC_vh_dt_phi.pdf : notebooks/lipid_xrr.py src/models/mol_vol.py
	cd notebooks && ipython lipid_xrr.py dppc 15 15 20 25 30 c
$(DMPG_OUT) $(FIG_DIR)/DMPG_ref_sld.pdf $(FIG_DIR)/DMPG_vh_dt_phi.pdf : notebooks/lipid_xrr.py src/models/mol_vol.py
	cd notebooks && ipython lipid_xrr.py dmpg 13 15 20 25 30 d
$(NDMPC1_OUT) : notebooks/lipid_nr.py $(DMPC_OUT)
	cd notebooks && ipython lipid_nr.py dmpc 13 25
$(NDMPC2_OUT) : notebooks/lipid_nr.py $(DMPC_OUT)
	cd notebooks && ipython lipid_nr.py dmpc 15 20
$(NDPPC1_OUT) : notebooks/lipid_nr.py $(DPPC_OUT)
	cd notebooks && ipython lipid_nr.py dppc 15 20
$(NDPPC2_OUT) : notebooks/lipid_nr.py $(DPPC_OUT)
	cd notebooks && ipython lipid_nr.py dppc 15 15


notebooks/lipid_xrr.py : notebooks/lipid_xrr.ipynb
	jupyter-nbconvert notebooks/lipid_xrr.ipynb --to script --output-dir=notebooks/
notebooks/lipid_nr.py : notebooks/lipid_nr.ipynb
	jupyter-nbconvert notebooks/lipid_nr.ipynb --to script --output-dir=notebooks/


$(FIG_DIR)/dmpc/20_neutron_corner.pdf : $(PLOT_DIR)/plotcornern.py output/dmpc20_chain_neutron.txt
	cd src/visualization && ipython plotcornern.py ../../output/dmpc20_chain_neutron.txt 13 dmpc 20
$(FIG_DIR)/dmpc/25_neutron_corner.pdf : $(PLOT_DIR)/plotcornern.py output/dmpc25_chain_neutron.txt
	cd src/visualization && ipython plotcornern.py ../../output/dmpc25_chain_neutron.txt 13 dmpc 25
$(FIG_DIR)/dppc/20_neutron_corner.pdf : $(PLOT_DIR)/plotcornern.py output/dppc20_chain_neutron.txt
	cd src/visualization && ipython plotcornern.py ../../output/dppc20_chain_neutron.txt 15 dppc 20
$(FIG_DIR)/dppc/15_neutron_corner.pdf : $(PLOT_DIR)/plotcornern.py output/dppc15_chain_neutron.txt
	cd src/visualization && ipython plotcornern.py ../../output/dppc15_chain_neutron.txt 15 dppc 15

$(PLOT_DIR)/plotcornern.py : $(PLOT_DIR)/plotcornern.ipynb
	jupyter-nbconvert src/visualization/plotcornern.ipynb --to script --output-dir=src/visualization/
