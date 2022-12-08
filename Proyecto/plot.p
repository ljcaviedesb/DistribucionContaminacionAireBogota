//No se como configurar el archivo para que salga bien la automatización del plot :/
set grid
set output 'Advec_Dif.png'
set title 'Plot de Adv_Dif'
splot "Ondas.dat" w l
