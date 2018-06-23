reports/paper.pdf : reports/paper.tex reports/rsc.bib reports/esi.pdf
	cd reports && pdflatex paper.tex
	cd reports && bibtex paper.aux
	cd reports && pdflatex paper.tex
reports/esi.pdf : reports/esi.tex
	cd reports && pdflatex esi.tex
	cd reports && bibtex esi.aux
	cd reports && pdflatex esi.tex
reports/paper.tex : reports/figures/* 
reports/figures/DLPC_all_data.png : output/dlpc4_sld.txt output/dlpc4_sld.txt output/dlpc4_sld.txt output/dlpc4_sld.txt output/dlpc_highconc_chain.txt src/visualization/plotref.py
	cd src/visualization && ipython plotref.py ../../output/dlpc4_ref.txt ../../output/dlpc4_sld.txt ../../output/dlpc5_ref.txt ../../output/dlpc5_sld.txt ../../output/dlpc_highconc_chain.txt a DLPC
reports/figures/dlpc4_all_corner.png : src/visualization/plotcorner.py output/dlpc_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dlpc_highconc_chain.txt 15.455 667 dlpc
output/dlpc_highconc_chain.txt : notebooks/DLPC/DLPC_all_conc.py
	cd notebooks/DLPC && ipython DLPC_all_conc.py
notebooks/DLPC/DLPC_all_conc.py : notebooks/DLPC/DLPC_all_conc.ipynb
	jupyter-nbconvert notebooks/DLPC/DLPC_all_conc.ipynb --to script --output-dir=notebooks/DLPC/
reports/figures/DMPC_all_data.png : src/visualization/plotref.py output/dmpc4_sld.txt output/dmpc4_sld.txt output/dmpc4_sld.txt output/dmpc4_sld.txt output/dmpc_highconc_chain.txt
	cd src/visualization && ipython plotref.py ../../output/dmpc4_ref.txt ../../output/dmpc4_sld.txt ../../output/dmpc5_ref.txt ../../output/dmpc5_sld.txt ../../output/dmpc_highconc_chain.txt b DMPC
reports/figures/dmpc4_all_corner.png : src/visualization/plotcorner.py output/dmpc_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dmpc_highconc_chain.txt 17.985 779 dmpc
output/dmpc_highconc_chain.txt : notebooks/DMPC/DMPC_all_conc.py
	cd notebooks/DMPC && ipython DMPC_all_conc.py
notebooks/DMPC/DMPC_all_conc.py : notebooks/DMPC/DMPC_all_conc.ipynb
	jupyter-nbconvert notebooks/DMPC/DMPC_all_conc.ipynb --to script --output-dir=notebooks/DMPC/
reports/figures/DPPC_all_data.png : output/dppc4_sld.txt output/dppc4_sld.txt output/dppc4_sld.txt output/dppc4_sld.txt output/dppc_highconc_chain.txt
	cd src/visualization && ipython plotref.py ../../output/dppc4_ref.txt ../../output/dppc4_sld.txt ../../output/dppc5_ref.txt ../../output/dppc5_sld.txt ../../output/dppc_highconc_chain.txt c DPPC
reports/figures/dppc4_all_corner.png : src/visualization/plotcorner.py output/dppc_highconc_chain.txt
	cd src/visualization && ipython plotcorner.py ../../output/dppc_highconc_chain.txt 20.515 891 dppc
output/dppc_highconc_chain.txt : notebooks/DPPC/DPPC_all_conc.py
	cd notebooks/DPPC && ipython DPPC_all_conc.py
notebooks/DPPC/DPPC_all_conc.py : notebooks/DPPC/DPPC_all_conc.ipynb
	jupyter-nbconvert notebooks/DPPC/DPPC_all_conc.ipynb --to script --output-dir=notebooks/DPPC/
reports/figures/DMPG_all_data.png : src/visualization/plotref.py output/dmpg4_sld.txt output/dmpg4_sld.txt output/dmpg4_sld.txt output/dmpg4_sld.txt output/dmpg_highconc_chain.txt
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

