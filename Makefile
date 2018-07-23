current_dir = $(shell pwd)
FIG_DIR = reports/figures
PLOT_DIR = src/visualization
DLPC_DIR = output/dlpc
DMPC_DIR = output/dmpc
DPPC_DIR = output/dppc
DMPG_DIR = output/dmpg
DATA_DIR = data/processed

PAPER_FIGS = $(FIG_DIR)/DLPC_all_data.png $(FIG_DIR)/DMPC_all_data.png $(FIG_DIR)/DMPG_all_data.png $(FIG_DIR)/DPPC_all_data.png $(FIG_DIR)/nDPPC_all_data.png $(FIG_DIR)/nDMPC_all_data.png $(FIG_DIR)/head_groups.png
ESI_FIGS = $(FIG_DIR)/dlpc4_all_corner.png $(FIG_DIR)/dlpc5_all_corner.png $(FIG_DIR)/dmpc3_neutron_corner1.png $(FIG_DIR)/dmpc3_neutron_corner2.png $(FIG_DIR)/dmpc4_all_corner.png $(FIG_DIR)/dmpc5_all_corner.png $(FIG_DIR)/dmpg4_all_corner.png $(FIG_DIR)/dmpg5_all_corner.png $(FIG_DIR)/dppc3_neutron_corner1.png $(FIG_DIR)/dppc3_neutron_corner2.png $(FIG_DIR)/dppc4_all_corner.png $(FIG_DIR)/dppc5_all_corner.png
DLPC_OUT = output/dlpc4_ref.txt output/dlpc4_sld.txt output/dlpc5_ref.txt output/dlpc5_sld.txt output/dlpc_highconc_chain.txt $(DLPC_DIR)/angle4.txt $(DLPC_DIR)/angle5.txt $(DLPC_DIR)/head4.txt $(DLPC_DIR)/head5.txt $(DLPC_DIR)/rough4.txt $(DLPC_DIR)/rough5.txt $(DLPC_DIR)/scale4.txt $(DLPC_DIR)/scale5.txt $(DLPC_DIR)/solh4.txt $(DLPC_DIR)/solh5.txt $(DLPC_DIR)/tail4.txt $(DLPC_DIR)/tail5.txt $(DLPC_DIR)/vh.txt
DMPC_OUT = output/dmpc4_ref.txt output/dmpc4_sld.txt output/dmpc5_ref.txt output/dmpc5_sld.txt output/dmpc_highconc_chain.txt $(DMPC_DIR)/angle4.txt $(DMPC_DIR)/angle5.txt $(DMPC_DIR)/head4.txt $(DMPC_DIR)/head5.txt $(DMPC_DIR)/rough4.txt $(DMPC_DIR)/rough5.txt $(DMPC_DIR)/scale4.txt $(DMPC_DIR)/scale5.txt $(DMPC_DIR)/solh4.txt $(DMPC_DIR)/solh5.txt $(DMPC_DIR)/tail4.txt $(DMPC_DIR)/tail5.txt $(DMPC_DIR)/vh.txt
DPPC_OUT = output/dppc4_ref.txt output/dppc4_sld.txt output/dppc5_ref.txt output/dppc5_sld.txt output/dppc_highconc_chain.txt $(DPPC_DIR)/angle4.txt $(DPPC_DIR)/angle5.txt $(DPPC_DIR)/head4.txt $(DPPC_DIR)/head5.txt $(DPPC_DIR)/rough4.txt $(DPPC_DIR)/rough5.txt $(DPPC_DIR)/scale4.txt $(DPPC_DIR)/scale5.txt $(DPPC_DIR)/solh4.txt $(DPPC_DIR)/solh5.txt $(DPPC_DIR)/tail4.txt $(DPPC_DIR)/tail5.txt $(DPPC_DIR)/vh.txt
DMPG_OUT = output/dmpg4_ref.txt output/dmpg4_sld.txt output/dmpg5_ref.txt output/dmpg5_sld.txt output/dmpg_highconc_chain.txt $(DMPG_DIR)/angle4.txt $(DMPG_DIR)/angle5.txt $(DMPG_DIR)/head4.txt $(DMPG_DIR)/head5.txt $(DMPG_DIR)/rough4.txt $(DMPG_DIR)/rough5.txt $(DMPG_DIR)/scale4.txt $(DMPG_DIR)/scale5.txt $(DMPG_DIR)/solh4.txt $(DMPG_DIR)/solh5.txt $(DMPG_DIR)/tail4.txt $(DMPG_DIR)/tail5.txt $(DMPG_DIR)/vh.txt
NDMPC1_OUT = output/dmpc3_n1_ref_neutron.txt output/dmpc3_n1_sld_neutron.txt output/dmpc_highconc_chain_neutron_n1.txt $(DMPC_DIR)/angle3_neutron_n1.txt $(DMPC_DIR)/rought3_neutron_n1.txt $(DMPC_DIR)/scale3_neutron_n1.txt $(DMPC_DIR)/solh3_neutron_n1.txt $(DMPC_DIR)/tail3_neutron_n1.txt
NDMPC2_OUT = output/dmpc3_n2_ref_neutron.txt output/dmpc3_n2_sld_neutron.txt output/dmpc_highconc_chain_neutron_n2.txt $(DMPC_DIR)/angle3_neutron_n2.txt $(DMPC_DIR)/rought3_neutron_n2.txt $(DMPC_DIR)/scale3_neutron_n2.txt $(DMPC_DIR)/solh3_neutron_n2.txt $(DMPC_DIR)/tail3_neutron_n2.txt
NDPPC1_OUT = output/dppc3_n1_ref_neutron.txt output/dppc3_n1_sld_neutron.txt output/dppc_highconc_chain_neutron_n1.txt $(DPPC_DIR)/angle3_neutron_n1.txt $(DPPC_DIR)/rought3_neutron_n1.txt $(DPPC_DIR)/scale3_neutron_n1.txt $(DPPC_DIR)/solh3_neutron_n1.txt $(DPPC_DIR)/tail3_neutron_n1.txt
NDPPC2_OUT = output/dppc3_n2_ref_neutron.txt output/dppc3_n2_sld_neutron.txt output/dppc_highconc_chain_neutron_n2.txt $(DPPC_DIR)/angle3_neutron_n2.txt $(DPPC_DIR)/rought3_neutron_n2.txt $(DPPC_DIR)/scale3_neutron_n2.txt $(DPPC_DIR)/solh3_neutron_n2.txt $(DPPC_DIR)/tail3_neutron_n2.txt

all : reports/paper.pdf reports/esi.pdf
clean :
	rm reports/paper.pdf reports/esi.pdf $(FIG_DIR)/*all_data.png $(FIG_DIR)/*all_corner.png $(FIG_DIR)/*neutron_corner*.png $(PLOT_DIR)/*.py notebooks/DLPC/*.py notebooks/DMPC/*.py notebooks/DPPC/*.py notebooks/DMPG/*.py output/* output/dlpc/* output/dppc/* output/dmpc/* output/dmpg/*

reports/paper.pdf : reports/paper.tex reports/rsc.bib $(PAPER_FIGS) $(DLPC_OUT) $(DMPC_OUT) $(DPPC_OUT) $(DMPG_OUT)
	cd reports && pdflatex paper.tex
	cd reports && bibtex paper.aux
	cd reports && pdflatex paper.tex
reports/esi.pdf : reports/esi.tex $(ESI_FIGS)
	cd reports && pdflatex esi.tex
	cd reports && bibtex esi.aux
	cd reports && pdflatex esi.tex


$(FIG_DIR)/DLPC_all_data.png : $(PLOT_DIR)/plotref.py $(DLPC_OUT)
	cd src/visualization && ipython plotref.py ../../output/dlpc4_ref.txt ../../output/dlpc4_sld.txt ../../output/dlpc5_ref.txt ../../output/dlpc5_sld.txt ../../output/dlpc_highconc_chain.txt a DLPC
$(FIG_DIR)/DMPC_all_data.png : $(PLOT_DIR)/plotref.py $(DMPC_OUT)
	cd src/visualization && ipython plotref.py ../../output/dmpc4_ref.txt ../../output/dmpc4_sld.txt ../../output/dmpc5_ref.txt ../../output/dmpc5_sld.txt ../../output/dmpc_highconc_chain.txt b DMPC
$(FIG_DIR)/DPPC_all_data.png : $(PLOT_DIR)/plotref.py $(DPPC_OUT)
	cd src/visualization && ipython plotref.py ../../output/dppc4_ref.txt ../../output/dppc4_sld.txt ../../output/dppc5_ref.txt ../../output/dppc5_sld.txt ../../output/dppc_highconc_chain.txt c DPPC
$(FIG_DIR)/DMPG_all_data.png : $(PLOT_DIR)/plotref.py $(DMPG_OUT)
	cd src/visualization && ipython plotref.py ../../output/dmpg4_ref.txt ../../output/dmpg4_sld.txt ../../output/dmpg5_ref.txt ../../output/dmpg5_sld.txt ../../output/dmpg_highconc_chain.txt d DMPG
$(FIG_DIR)/nDMPC_all_data.png : $(PLOT_DIR)/plotrefn.py $(NDMPC1_OUT) $(NDMPC2_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dmpc3_n1_ref_neutron.txt ../../output/dmpc3_n1_sld_neutron.txt ../../output/dmpc3_n2_ref_neutron.txt ../../output/dmpc3_n2_sld_neutron.txt ../../output/dmpc_highconc_chain_neutron_n2.txt a nDMPC
$(FIG_DIR)/nDPPC_all_data.png : $(PLOT_DIR)/plotrefn.py $(NDPPC1_OUT) $(NDPPC2_OUT)
	cd src/visualization && ipython plotrefn.py ../../output/dppc3_n1_ref_neutron.txt ../../output/dppc3_n1_sld_neutron.txt ../../output/dppc3_n2_ref_neutron.txt ../../output/dppc3_n2_sld_neutron.txt ../../output/dppc_highconc_chain_neutron_n2.txt b nDPPC

$(PLOT_DIR)/plotref.py : $(PLOT_DIR)/plotref.ipynb
	jupyter-nbconvert src/visualization/plotref.ipynb --to script --output-dir=src/visualization/
$(PLOT_DIR)/plotrefn.py : $(PLOT_DIR)/plotrefn.ipynb
	jupyter-nbconvert src/visualization/plotrefn.ipynb --to script --output-dir=src/visualization/

$(DLPC_OUT) : notebooks/DLPC/DLPC_all_conc.py $(DATA_DIR)/DLPC/DLPC_Xray_conc4.dat $(DATA_DIR)/DLPC/DLPC_Xray_conc5.dat
	cd notebooks/DLPC && ipython DLPC_all_conc.py $(current_dir)
$(DMPC_OUT) : notebooks/DMPC/DMPC_all_conc.py $(DATA_DIR)/DMPC/DMPC_Xray_conc4.dat $(DATA_DIR)/DMPC/DMPC_Xray_conc5.dat
	cd notebooks/DMPC && ipython DMPC_all_conc.py $(current_dir)
$(DPPC_OUT) : notebooks/DPPC/DPPC_all_conc.py $(DATA_DIR)/DPPC/DPPC_Xray_conc4.dat $(DATA_DIR)/DPPC/DPPC_Xray_conc5.dat
	cd notebooks/DPPC && ipython DPPC_all_conc.py $(current_dir)
$(DMPG_OUT) : notebooks/DMPG/DMPG_all_conc.py $(DATA_DIR)/DMPG/DMPG_Xray_conc4.dat $(DATA_DIR)/DMPG/DMPG_Xray_conc5.dat
	cd notebooks/DMPG && ipython DMPG_all_conc.py $(current_dir)
$(NDMPC1_OUT) : notebooks/DMPC/DMPC_neutron_conc3_n1.py $(DATA_DIR)/DMPC/DMPC_Neutron_conc2_dDMPC_hDES.mft
	cd notebooks/DMPC && ipython DMPC_neutron_conc3_n1.py $(current_dir)
$(NDMPC2_OUT) : notebooks/DMPC/DMPC_neutron_conc3_n2.py $(DATA_DIR)/DMPC/DMPC_Neutron_conc1_dDMPC_hdDES.mft
	cd notebooks/DMPC && ipython DMPC_neutron_conc3_n2.py $(current_dir)
$(NDPPC1_OUT) : notebooks/DPPC/DPPC_neutron_conc3_n1.py $(DATA_DIR)/DPPC/DPPC_Neutron_conc3_dDPPC_hDES.mft
	cd notebooks/DPPC && ipython DPPC_neutron_conc3_n1.py $(current_dir)
$(NDPPC2_OUT) : notebooks/DPPC/DPPC_neutron_conc3_n2.py $(DATA_DIR)/DPPC/DPPC_Neutron_conc3_dDPPC_hdDES.mft
	cd notebooks/DPPC && ipython DPPC_neutron_conc3_n2.py $(current_dir)

notebooks/DLPC/DLPC_all_conc.py : notebooks/DLPC/DLPC_all_conc.ipynb
	jupyter-nbconvert notebooks/DLPC/DLPC_all_conc.ipynb --to script --output-dir=notebooks/DLPC/
notebooks/DMPC/DMPC_all_conc.py : notebooks/DMPC/DMPC_all_conc.ipynb
	jupyter-nbconvert notebooks/DMPC/DMPC_all_conc.ipynb --to script --output-dir=notebooks/DMPC/
notebooks/DPPC/DPPC_all_conc.py : notebooks/DPPC/DPPC_all_conc.ipynb
	jupyter-nbconvert notebooks/DPPC/DPPC_all_conc.ipynb --to script --output-dir=notebooks/DPPC/
notebooks/DMPG/DMPG_all_conc.py : notebooks/DMPG/DMPG_all_conc.ipynb
	jupyter-nbconvert notebooks/DMPG/DMPG_all_conc.ipynb --to script --output-dir=notebooks/DMPG/
notebooks/DMPC/DMPC_neutron_conc3_n1.py : notebooks/DMPC/DMPC_neutron_conc3_n1.ipynb
	jupyter-nbconvert notebooks/DMPC/DMPC_neutron_conc3_n1.ipynb --to script --output-dir=notebooks/DMPC/
notebooks/DMPC/DMPC_neutron_conc3_n2.py : notebooks/DMPC/DMPC_neutron_conc3_n2.ipynb
	jupyter-nbconvert notebooks/DMPC/DMPC_neutron_conc3_n2.ipynb --to script --output-dir=notebooks/DMPC/
notebooks/DPPC/DPPC_neutron_conc3_n1.py : notebooks/DPPC/DPPC_neutron_conc3_n1.ipynb
	jupyter-nbconvert notebooks/DPPC/DPPC_neutron_conc3_n1.ipynb --to script --output-dir=notebooks/DPPC/
notebooks/DPPC/DPPC_neutron_conc3_n2.py : notebooks/DPPC/DPPC_neutron_conc3_n2.ipynb
	jupyter-nbconvert notebooks/DPPC/DPPC_neutron_conc3_n2.ipynb --to script --output-dir=notebooks/DPPC/



$(FIG_DIR)/dlpc4_all_corner.png $(FIG_DIR)/dlpc5_all_corner.png : $(PLOT_DIR)/plotcorner.py output/dlpc_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dlpc_highconc_chain.txt 15.455 667 dlpc $(current_dir)
$(FIG_DIR)/dmpc4_all_corner.png $(FIG_DIR)/dmpc5_all_corner.png : $(PLOT_DIR)/plotcorner.py output/dmpc_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dmpc_highconc_chain.txt 17.985 779 dmpc $(current_dir)
$(FIG_DIR)/dppc4_all_corner.png $(FIG_DIR)/dppc5_all_corner.png : $(PLOT_DIR)/plotcorner.py output/dppc_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dppc_highconc_chain.txt 20.515 891 dppc $(current_dir)
$(FIG_DIR)/dmpg4_all_corner.png $(FIG_DIR)/dmpg5_all_corner.png : $(PLOT_DIR)/plotcorner.py output/dmpg_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dmpg_highconc_chain.txt 17.985 779 dmpg $(current_dir)


$(FIG_DIR)/dmpc3_neutron_corner1.png : $(PLOT_DIR)/plotcornern.py output/dmpc_highconc_chain_neutron_n1.txt
	cd src/visualization && ipython plotcornern.py ../../output/dmpc_highconc_chain_neutron_n1.txt 17.985 779 dmpc 1 $(current_dir)
$(FIG_DIR)/dmpc3_neutron_corner2.png : $(PLOT_DIR)/plotcornern.py output/dmpc_highconc_chain_neutron_n2.txt
	cd src/visualization && ipython plotcornern.py ../../output/dmpc_highconc_chain_neutron_n2.txt 17.985 779 dmpc 2 $(current_dir)
$(FIG_DIR)/dppc3_neutron_corner1.png : $(PLOT_DIR)/plotcornern.py output/dppc_highconc_chain_neutron_n1.txt
	cd src/visualization && ipython plotcornern.py ../../output/dppc_highconc_chain_neutron_n1.txt 20.515 891 dppc 1 $(current_dir)
$(FIG_DIR)/dppc3_neutron_corner2.png : $(PLOT_DIR)/plotcornern.py output/dppc_highconc_chain_neutron_n2.txt
	cd src/visualization && ipython plotcornern.py ../../output/dppc_highconc_chain_neutron_n2.txt 20.515 891 dppc 2 $(current_dir)


$(PLOT_DIR)/plotcorner.py : $(PLOT_DIR)/plotcorner.ipynb
	jupyter-nbconvert src/visualization/plotcorner.ipynb --to script --output-dir=src/visualization/
$(PLOT_DIR)/plotcornern.py : $(PLOT_DIR)/plotcornern.ipynb
	jupyter-nbconvert src/visualization/plotcornern.ipynb --to script --output-dir=src/visualization/

