figure-neuroblastomaProcessed-combinations.png: figure-neuroblastomaProcessed-combinations.R
	R --vanilla < $<
figure-neuroblastomaProcessed-complex/index.html: figure-neuroblastomaProcessed-complex.R
	R --vanilla < $<

