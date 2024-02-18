pre: runmain plot

runmain:
	python PF_main_onecomponent.py
	python PF_main_multicomponent.py
	
plot:
	python plot_histogram_wmax.py

