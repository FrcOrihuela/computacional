#
# Representacion Tcuad vs L
# Practica pendulo simple
#

#siempre estará el dibujo en la pantalla para poder hacer otra cosa
set terminal
#Título de los ejes
set xlabel "Voltaje (V)})"
set ylabel "Intensidad (nA)" enhanced

#dibujo sin leyenda
set nokey

#y(x) será la función a ajustar
a=1.0
b=1.0
c=1.0
y(x)=b*exp(-a*x)

# ajusta los datos del fichero mediante los parámetros a y b
fit y(x) "barreras.txt" u 1:2:3 yerrors via a, b 
#dibuja las columnas del archivo "datos.txt" y la función a ajustar 
plot  "barreras.txt" u 1:2:3 w yerrors, y(x)

set terminal postscript eps enhanced color 14

#en vez de mandar a dibujar a windows que no lo entiende y muestra todos los comandos, esta instrucción hace que se almacene en un fichero:
set output "barreras.png"

#para hacerlo uniforme sobre todo el papel, aprovechar toda la escala
set size 1,1

#redibujar
replot


set terminal