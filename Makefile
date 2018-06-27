reports/paper.pdf : reports/paper.tex reports/rsc.bib reports/esi.pdf reports/figures/* output/* output/dlpc/* output/dppc/* output/dmpc/* output/dmpg/*
	cd reports && pdflatex paper.tex
	cd reports && bibtex paper.aux
	cd reports && pdflatex paper.tex
reports/esi.pdf : reports/esi.tex reports/figures/*
	cd reports && pdflatex esi.tex
	cd reports && bibtex esi.aux
	cd reports && pdflatex esi.tex
reports/figures/DLPC_all_data.png : output/dlpc4_sld.txt output/dlpc4_ref.txt output/dlpc5_sld.txt output/dlpc5_ref.txt output/dlpc_highconc_chain.txt src/visualization/plotref.py
	cd src/visualization && ipython plotref.py ../../output/dlpc4_ref.txt ../../output/dlpc4_sld.txt ../../output/dlpc5_ref.txt ../../output/dlpc5_sld.txt ../../output/dlpc_highconc_chain.txt a DLPC
reports/figures/dlpc4_all_corner.png : src/visualization/plotcorner.py output/dlpc_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dlpc_highconc_chain.txt 15.455 667 dlpc
output/dlpc_highconc_chain.txt : notebooks/DLPC/DLPC_all_conc.py
	cd notebooks/DLPC && ipython DLPC_all_conc.py
notebooks/DLPC/DLPC_all_conc.py : notebooks/DLPC/DLPC_all_conc.ipynb
	jupyter-nbconvert notebooks/DLPC/DLPC_all_conc.ipynb --to script --output-dir=notebooks/DLPC/
reports/figures/DMPC_all_data.png : src/visualization/plotref.py output/dmpc4_ref.txt output/dmpc4_sld.txt output/dmpc5_ref.txt output/dmpc5_sld.txt output/dmpc_highconc_chain.txt
	cd src/visualization && ipython plotref.py ../../output/dmpc4_ref.txt ../../output/dmpc4_sld.txt ../../output/dmpc5_ref.txt ../../output/dmpc5_sld.txt ../../output/dmpc_highconc_chain.txt b DMPC
reports/figures/nDMPC_all_data.png : output/dmpc3_n1_sld_neutron.txt output/dmpc3_n1_ref_neutron.txt output/dmpc3_n2_sld_neutron.txt output/dmpc3_n2_ref_neutron.txt 
	cd src/visualization && ipython plotrefn.py ../../output/dmpc3_n1_ref_neutron.txt ../../output/dmpc3_n1_sld_neutron.txt ../../output/dmpc3_n2_ref_neutron.txt ../../output/dmpc3_n2_sld_neutron.txt ../../output/dmpc_highconc_chain_neutron_n2.txt a nDMPC
reports/figures/dmpc4_all_corner.png : src/visualization/plotcorner.py output/dmpc_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dmpc_highconc_chain.txt 17.985 779 dmpc
reports/figures/dmpc3_neutron_corner1.png : src/visualization/plotcornern.py output/dmpc_highconc_chain_neutron_n1.txt
	cd src/visualization && ipython plotcornern.py ../../output/dmpc_highconc_chain_neutron_n1.txt 17.985 779 dmpc 1
reports/figures/dmpc3_neutron_corner2.png : src/visualization/plotcornern.py output/dmpc_highconc_chain_neutron_n2.txt
	cd src/visualization && ipython plotcornern.py ../../output/dmpc_highconc_chain_neutron_n2.txt 17.985 779 dmpc 2
output/dmpc_highconc_chain.txt : notebooks/DMPC/DMPC_all_conc.py
	cd notebooks/DMPC && ipython DMPC_all_conc.py
notebooks/DMPC/DMPC_all_conc.py : notebooks/DMPC/DMPC_all_conc.ipynb
	jupyter-nbconvert notebooks/DMPC/DMPC_all_conc.ipynb --to script --output-dir=notebooks/DMPC/
output/dmpc3_highconc_chain_neutron_n1.txt : notebooks/DMPC/DMPC_neutron_conc3_n1.py
	cd notebooks/DMPC && ipython DMPC_neutron_conc3_n1.py
notebooks/DMPC/DMPC_neutron_conc3_n1.py : notebooks/DMPC/DMPC_neutron_conc3_n1.ipynb
	jupyter-nbconvert notebooks/DMPC/DMPC_neutron_conc3_n1.ipynb --to script --output-dir=notebooks/DMPC/
output/dmpc3_highconc_chain_neutron_n2.txt : notebooks/DMPC/DMPC_neutron_conc3_n2.py
	cd notebooks/DMPC && ipython DMPC_neutron_conc3_n2.py
notebooks/DMPC/DMPC_neutron_conc3_n2.py : notebooks/DMPC/DMPC_neutron_conc3_n2.ipynb
	jupyter-nbconvert notebooks/DMPC/DMPC_neutron_conc3_n2.ipynb --to script --output-dir=notebooks/DMPC/
reports/figures/DPPC_all_data.png : output/dppc4_sld.txt output/dppc4_ref.txt output/dppc5_sld.txt output/dppc5_ref.txt output/dppc_highconc_chain.txt
	cd src/visualization && ipython plotref.py ../../output/dppc4_ref.txt ../../output/dppc4_sld.txt ../../output/dppc5_ref.txt ../../output/dppc5_sld.txt ../../output/dppc_highconc_chain.txt c DPPC
reports/figures/nDPPC_all_data.png : output/dppc3_n1_sld_neutron.txt output/dppc3_n1_ref_neutron.txt output/dppc3_n2_sld_neutron.txt output/dppc3_n2_ref_neutron.txt 
	cd src/visualization && ipython plotrefn.py ../../output/dppc3_n1_ref_neutron.txt ../../output/dppc3_n1_sld_neutron.txt ../../output/dppc3_n2_ref_neutron.txt ../../output/dppc3_n2_sld_neutron.txt ../../output/dppc_highconc_chain_neutron_n2.txt b nDPPC
reports/figures/dppc4_all_corner.png : src/visualization/plotcorner.py output/dppc_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dppc_highconc_chain.txt 20.515 891 dppc
reports/figures/dppc3_neutron_corner1.png : src/visualization/plotcornern.py output/dppc_highconc_chain_neutron_n1.txt
	cd src/visualization && ipython plotcornern.py ../../output/dppc_highconc_chain_neutron_n1.txt 20.515 891 dppc 1
reports/figures/dppc3_neutron_corner2.png : src/visualization/plotcornern.py output/dppc_highconc_chain_neutron_n2.txt
	cd src/visualization && ipython plotcornern.py ../../output/dppc_highconc_chain_neutron_n2.txt 20.515 891 dppc 2
output/dppc_highconc_chain.txt : notebooks/DPPC/DPPC_all_conc.py
	cd notebooks/DPPC && ipython DPPC_all_conc.py
notebooks/DPPC/DPPC_all_conc.py : notebooks/DPPC/DPPC_all_conc.ipynb
	jupyter-nbconvert notebooks/DPPC/DPPC_all_conc.ipynb --to script --output-dir=notebooks/DPPC/
output/dppc3_highconc_chain_neutron_n1.txt : notebooks/DPPC/DPPC_neutron_conc3_n1.py
	cd notebooks/DPPC && ipython DPPC_neutron_conc3_n1.py
notebooks/DPPC/DPPC_neutron_conc3_n1.py : notebooks/DPPC/DPPC_neutron_conc3_n1.ipynb
	jupyter-nbconvert notebooks/DPPC/DPPC_neutron_conc3_n1.ipynb --to script --output-dir=notebooks/DPPC/
output/dppc3_highconc_chain_neutron_n2.txt : notebooks/DPPC/DPPC_neutron_conc3_n2.py
	cd notebooks/DPPC && ipython DPPC_neutron_conc3_n2.py
notebooks/DPPC/DPPC_neutron_conc3_n2.py : notebooks/DPPC/DPPC_neutron_conc3_n2.ipynb
	jupyter-nbconvert notebooks/DPPC/DPPC_neutron_conc3_n2.ipynb --to script --output-dir=notebooks/DPPC/
reports/figures/dppc3_neutron_corner.png : src/visualization/plotcornern.py output/dppc_highconc_chain_neutron.txt
	cd src/visualization && ipython plotcornern.py ../../output/dppc_highconc_chain_neutron.txt 20.515 891 dppc
reports/figures/DMPG_all_data.png : src/visualization/plotref.py output/dmpg4_sld.txt output/dmpg4_ref.txt output/dmpg5_sld.txt output/dmpg5_ref.txt output/dmpg_highconc_chain.txt
	cd src/visualization && ipython plotref.py ../../output/dmpg4_ref.txt ../../output/dmpg4_sld.txt ../../output/dmpg5_ref.txt ../../output/dmpg5_sld.txt ../../output/dmpg_highconc_chain.txt d DMPG
reports/figures/dmpg4_all_corner.png : src/visualization/plotcorner.py output/dmpg_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dmpg_highconc_chain.txt 17.985 779 dmpg
output/dmpg_highconc_chain.txt : notebooks/DMPG/DMPG_all_conc.py
	cd notebooks/DMPG && ipython DMPG_all_conc.py
notebooks/DMPG/DMPG_all_conc.py : notebooks/DMPG/DMPG_all_conc.ipynb
	jupyter-nbconvert notebooks/DMPG/DMPG_all_conc.ipynb --to script --output-dir=notebooks/DMPG/
src/visualization/plotref.py : src/visualization/plotref.ipynb
	jupyter-nbconvert src/visualization/plotref.ipynb --to script -- --output-dir=src/visualization
src/visualization/plotcorner.py : src/visualization/plotcorner.ipynb
	jupyter-nbconvert src/visualization/plotcorner.ipynb --to script -- --output-dir=src/visualization
src/visualization/plotrefn.py : src/visualization/plotrefn.ipynb
	jupyter-nbconvert src/visualization/plotrefn.ipynb --to script -- --output-dir=src/visualization
src/visualization/plotcornern.py : src/visualization/plotcornern.ipynb
	jupyter-nbconvert src/visualization/plotcornern.ipynb --to script -- --output-dir=src/visualization

