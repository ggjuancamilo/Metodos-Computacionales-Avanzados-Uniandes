evst.pdf : data.txt graficar.py
	python3 graficar.py

data.txt tiempo_vs_procesadores.pdf : tarea1.out Tiempo.py
	python3 Tiempo.py

tarea1.out : tarea1.c
	gcc -fopenmp tarea1.c -lm -o tarea1.out

clean:
	rm -f tarea1.out data.txt borrar.txt tiempo1.txt tiempo2.txt tiempo3.txt
