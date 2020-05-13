exec:
	gcc shock.c -lm -o exec_shock.x
	gcc funciones.c main.c -lm -o exec_sedov.x

shock:
	./exec_shock.x

sedov:
	./exec_sedov.x

plotshock:
	python3 plotter.py

plotsedov:
	python3 explot.py

clean:
	rm -f exec_shock.x exec_sedov.x datos data1 data2 data3
