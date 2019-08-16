figure-neuroblastomaProcessed-combinations-interactive/index.html: figure-neuroblastomaProcessed-combinations-interactive.R neuroblastomaProcessed.combinations.rds
	R --vanilla < $<
figure-neuroblastomaProcessed-combinations.png: figure-neuroblastomaProcessed-combinations.R neuroblastomaProcessed.combinations.rds
	R --vanilla < $<
figure-curveAlignment/index.html: figure-curveAlignment.R curveAlignment.rds
	R --vanilla < $<
curveAlignment.rds: curveAlignment.R
	R --vanilla < $<
neuroblastomaProcessed.combinations.rds: neuroblastomaProcessed.combinations.R
	R --vanilla < $<
figure-neuroblastomaProcessed-complex/index.html: figure-neuroblastomaProcessed-complex.R
	R --vanilla < $<

