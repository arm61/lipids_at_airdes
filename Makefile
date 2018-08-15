current_dir = $(shell pwd)
FIG_DIR = reports/figures
PLOT_DIR = src/visualization
DLPC_DIR = output/dlpc
DMPC_DIR = output/dmpc
DPPC_DIR = output/dppc
DMPG_DIR = output/dmpg
DATA_DIR = data/processed

PAPER_FIGS = $(FIG_DIR)/DLPC_all_data.pdf $(FIG_DIR)/DMPC_all_data.pdf $(FIG_DIR)/DMPG_all_data.pdf $(FIG_DIR)/DPPC_all_data.pdf $(FIG_DIR)/nDPPC15_all_data.pdf $(FIG_DIR)/nDMPC20_all_data.pdf $(FIG_DIR)/nDPPC20_all_data.pdf $(FIG_DIR)/nDMPC25_all_data.pdf $(FIG_DIR)/head_groups.png
ESI_FIGS = $(FIG_DIR)/dlpc4_all_corner.pdf $(FIG_DIR)/dlpc5_all_corner.pdf $(FIG_DIR)/dmpc25_neutron_corner.pdf $(FIG_DIR)/dmpc20_neutron_corner.pdf $(FIG_DIR)/dmpc4_all_corner.pdf $(FIG_DIR)/dmpc5_all_corner.pdf $(FIG_DIR)/dmpg4_all_corner.pdf $(FIG_DIR)/dmpg5_all_corner.pdf $(FIG_DIR)/dppc20_neutron_corner.pdf $(FIG_DIR)/dppc15_neutron_corner.pdf $(FIG_DIR)/dppc4_all_corner.pdf $(FIG_DIR)/dppc5_all_corner.pdf
DLPC_OUT = output/dlpc_chain.txt
DMPC_OUT = output/dmpc_chain.txt
DPPC_OUT = output/dppc_chain.txt
DMPG_OUT = output/dmpg_chain.txt
NDMPC1_OUT = output/dmpc25_chain_neutron.txt
NDMPC2_OUT = output/dmpc20_chain_neutron.txt
NDPPC1_OUT = output/dppc20_chain_neutron.txt
NDPPC2_OUT = output/dppc15_chain_neutron.txt

all : reports/paper.pdf reports/esi.pdf
clean :
	rm reports/paper.pdf reports/esi.pdf $(FIG_DIR)/*all_data.pdf $(FIG_DIR)/*all_corner.pdf $(FIG_DIR)/*neutron_corner*.pdf $(PLOT_DIR)/*.py notebooks/DLPC/*.py notebooks/DMPC/*.py notebooks/DPPC/*.py notebooks/DMPG/*.py output/* output/dlpc/* output/dppc/* output/dmpc/* output/dmpg/*

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


$(FIG_DIR)/DLPC_all_data.pdf : $(PLOT_DIR)/plotref.py $(DLPC_OUT)
	cd src/visualization && ipython plotref.py ../../output/dlpc20_ref.txt ../../output/dlpc20_sld.txt ../../output/dlpc25_ref.txt ../../output/dlpc25_sld.txt ../../output/dlpc30_ref.txt ../../output/dlpc30_sld.txt ../../output/dlpc35_ref.txt ../../output/dlpc35_sld.txt ../../output/dlpc_chain.txt a DLPC
$(FIG_DIR)/DMPC_all_data.pdf : $(PLOT_DIR)/plotref.py $(DMPC_OUT)
	cd src/visualization && ipython plotref.py ../../output/dmpc20_ref.txt ../../output/dmpc20_sld.txt ../../output/dmpc25_ref.txt ../../output/dmpc25_sld.txt ../../output/dmpc30_ref.txt ../../output/dmpc30_sld.txt ../../output/dmpc40_ref.txt ../../output/dmpc40_sld.txt ../../output/dmpc_chain.txt b DMPC
$(FIG_DIR)/DPPC_all_data.pdf : $(PLOT_DIR)/plotref.py $(DPPC_OUT)
	cd src/visualization && ipython plotref.py ../../output/dppc15_ref.txt ../../output/dppc15_sld.txt ../../output/dppc20_ref.txt ../../output/dppc20_sld.txt ../../output/dppc25_ref.txt ../../output/dppc25_sld.txt ../../output/dppc30_ref.txt ../../output/dppc30_sld.txt ../../output/dppc_chain.txt c DPPC
$(FIG_DIR)/DMPG_all_data.pdf : $(PLOT_DIR)/plotref.py $(DMPG_OUT)
	cd src/visualization && ipython plotref.py ../../output/dmpg15_ref.txt ../../output/dmpg15_sld.txt ../../output/dmpg20_ref.txt ../../output/dmpg20_sld.txt ../../output/dmpg25_ref.txt ../../output/dmpg25_sld.txt ../../output/dmpg30_ref.txt ../../output/dmpg30_sld.txt ../../output/dmpg_chain.txt d DMPG
$(FIG_DIR)/nDMPC25_all_data.pdf : $(PLOT_DIR)/plotrefn.py $(NDMPC1_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dmpc25_1_ref_neutron.txt ../../output/dmpc25_1_sld_neutron.txt ../../output/dmpc25_2_ref_neutron.txt ../../output/dmpc25_2_sld_neutron.txt ../../output/dmpc25_chain_neutron.txt a nDMPC25
$(FIG_DIR)/nDMPC20_all_data.pdf : $(PLOT_DIR)/plotrefn.py $(NDMPC2_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dmpc20_1_ref_neutron.txt ../../output/dmpc20_1_sld_neutron.txt ../../output/dmpc20_2_ref_neutron.txt ../../output/dmpc20_2_sld_neutron.txt ../../output/dmpc20_chain_neutron.txt a nDMPC20
$(FIG_DIR)/nDPPC20_all_data.pdf : $(PLOT_DIR)/plotrefn.py $(NDPPC1_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dppc20_1_ref_neutron.txt ../../output/dppc20_1_sld_neutron.txt ../../output/dppc20_2_ref_neutron.txt ../../output/dppc20_2_sld_neutron.txt ../../output/dppc20_chain_neutron.txt b nDPPC20
$(FIG_DIR)/nDPPC15_all_data.pdf : $(PLOT_DIR)/plotrefn.py $(NDPPC2_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dppc15_1_ref_neutron.txt ../../output/dppc15_1_sld_neutron.txt ../../output/dppc15_2_ref_neutron.txt ../../output/dppc15_2_sld_neutron.txt ../../output/dppc15_chain_neutron.txt b nDPPC15

$(PLOT_DIR)/plotref.py : $(PLOT_DIR)/plotref.ipynb
	jupyter-nbconvert src/visualization/plotref.ipynb --to script --output-dir=src/visualization/
$(PLOT_DIR)/plotrefn.py : $(PLOT_DIR)/plotrefn.ipynb
	jupyter-nbconvert src/visualization/plotrefn.ipynb --to script --output-dir=src/visualization/

$(DLPC_OUT) : notebooks/lipid_xrr.py src/models/mol_vol.py
	cd notebooks && ipython lipid_xrr.py dlpc 11 20 25 30 35
$(DMPC_OUT) : notebooks/lipid_xrr.py src/models/mol_vol.py
	cd notebooks && ipython lipid_xrr.py dmpc 13 20 25 30 40
$(DPPC_OUT) : notebooks/lipid_xrr.py src/models/mol_vol.py
	cd notebooks && ipython lipid_xrr.py dppc 15 15 20 25 30
$(DMPG_OUT) : notebooks/lipid_xrr.py src/models/mol_vol.py
	cd notebooks && ipython lipid_xrr.py dmpg 13 15 20 25 30
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


$(FIG_DIR)/dlpc4_all_corner.pdf $(FIG_DIR)/dlpc5_all_corner.pdf : $(PLOT_DIR)/plotcorner.py output/dlpc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dlpc_chain.txt 11 dlpc
$(FIG_DIR)/dmpc4_all_corner.pdf $(FIG_DIR)/dmpc5_all_corner.pdf : $(PLOT_DIR)/plotcorner.py output/dmpc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dmpc_chain.txt 13 dmpc
$(FIG_DIR)/dppc4_all_corner.pdf $(FIG_DIR)/dppc5_all_corner.pdf : $(PLOT_DIR)/plotcorner.py output/dppc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dppc_chain.txt 15 dppc
$(FIG_DIR)/dmpg4_all_corner.pdf $(FIG_DIR)/dmpg5_all_corner.pdf : $(PLOT_DIR)/plotcorner.py output/dmpg_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dmpg_chain.txt 13 dmpg


$(FIG_DIR)/dmpc20_neutron_corner.pdf : $(PLOT_DIR)/plotcornern.py output/dmpc20_chain_neutron.txt
	cd src/visualization && ipython plotcornern.py ../../output/dmpc20_chain_neutron.txt 13 dmpc 20
$(FIG_DIR)/dmpc25_neutron_corner.pdf : $(PLOT_DIR)/plotcornern.py output/dmpc25_chain_neutron.txt
	cd src/visualization && ipython plotcornern.py ../../output/dmpc25_chain_neutron.txt 13 dmpc 25
$(FIG_DIR)/dppc20_neutron_corner.pdf : $(PLOT_DIR)/plotcornern.py output/dppc20_chain_neutron.txt
	cd src/visualization && ipython plotcornern.py ../../output/dppc20_chain_neutron.txt 15 dppc 20
$(FIG_DIR)/dppc15_neutron_corner.pdf : $(PLOT_DIR)/plotcornern.py output/dppc15_chain_neutron.txt
	cd src/visualization && ipython plotcornern.py ../../output/dppc15_chain_neutron.txt 15 dppc 15


$(PLOT_DIR)/plotcorner.py : $(PLOT_DIR)/plotcorner.ipynb
	jupyter-nbconvert src/visualization/plotcorner.ipynb --to script --output-dir=src/visualization/
$(PLOT_DIR)/plotcorner-new.py : $(PLOT_DIR)/plotcorner-new.ipynb
	jupyter-nbconvert src/visualization/plotcorner-new.ipynb --to script --output-dir=src/visualization/
$(PLOT_DIR)/plotcornern.py : $(PLOT_DIR)/plotcornern.ipynb
	jupyter-nbconvert src/visualization/plotcornern.ipynb --to script --output-dir=src/visualization/
