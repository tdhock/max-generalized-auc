figure-curveAlignment/index.html: figure-curveAlignment.R curveAlignment.rds
	R --vanilla < $<
curveAlignment.rds: curveAlignment.R
	R --vanilla < $<
figure-neuroblastomaProcessed-combinations.png: figure-neuroblastomaProcessed-combinations.R
	R --vanilla < $<
figure-neuroblastomaProcessed-complex/index.html: figure-neuroblastomaProcessed-complex.R
	R --vanilla < $<

