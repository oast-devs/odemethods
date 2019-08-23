#include "ODE_methods.h"


/* 
 * ODE1 integrator is based on Euler integration (https://en.wikipedia.org/wiki/Euler_method)
 * the function takes the following inputs:
 *		x0: starting state of the differential equation
 *		t: end time of the integration
 *		stepsize: fixed size of the euler integration step
 *		output: pointer to the output vector
 *		fnc_ptr: pointer to the differential equation in the form of dx/dt = function(t, x)
 *
 * output:
 * 		int nsteps: number of integration steps performed  
 */

/* Vi spiego alcuni dettagli implementativi per principianti di C:
 * il primo dettaglio che vorrei evidenziare è il double **output.
 * si tratta di un puntatore ad un puntatore. Un puntatore ad un puntatore
 * può essere visto come il puntatore al puntatore di un'area di memoria, quindi come 
 * l'indirizzo di un vettore di elementi. Il doppio puntatore serve perché noi non possiamo
 * conoscere da principio la dimensione dell'output del problema (che invece ritorniamo
 * come valore di ritorno della funzione) ma riallochiamo la dimensione del vettore
 * all'interno della funzione stessa.
 * Per far si che possiamo correttamente riallocare il vettore e rendere visibile la
 * modifica al di fuori della funzione ODE, abbiamo bisogno di un doppio puntatore
 * che conterrà quindi l'indirizzo al nuovo vettore di output creato durante l'integrazione.
 *
 * L'altro dettaglio interessante è l'utilizzo di un puntatore a funzione.
 * Questo puntatore funziona similmente agli handle di matlab. In questo modo
 * la ODE rimane generica e la funzione differenziale da integrare viene passata
 * sotto forma di puntatore. La differenza fondamentale rispetto a Matlab è che in
 * C va dichiarata la forma della funzione (double funzione(double, double) nel
 * nostro caso) altrimenti il codice non compila.
 *
 * Futuri aggiornamenti al codice prevedono che la funzione da integrare sia vettoriale7
 * (cioè fnc_ptr può prendere in ingresso 3 stati e tornate 3 derivate x', y', z'). Per fare questo
 * le modifiche necessarie sono:
 *
 * double *x0 (vettore di stati iniziali)
 * double *** output (vettore di vettori di stati nel tempo. Con un opportuno protocollo si
 *                    potrebbe mantenere solo il doppio vettore in realtà)
 * double *(*fnc_ptr)(double, double *) (funzione da integrare che prende in ingresso il tempo
 *										 t e il vettore di stati e ritorno il vettore di derivate).
 *
 * Di seguito alcuni commenti illustreranno meglio il codice
 *
 */
int ODE1(double x0, double t, double stepsize, double **output, double(*fnc_ptr)(double, double))
{
	// controllo che tutti gli ingressi abbiano valori accettabili
	if(t < 0){
		fprintf(stderr, "Error, invalid time t\n");
		return -1;
	}
	
	if(stepsize <= 0) {
		fprintf(stderr, "Error, step size should be > 0\n");
		return -1;
	}
	
	if(output == NULL){
		fprintf(stderr, "Error, invalid output vector\n");
		return -1;
	}
	
	if(fnc_ptr == NULL) {
		fprintf(stderr, "Error, invalid ODE function\n");
		return -1;
	}
	
	// Calcolo il numero di step del problema, verificando
	// che sia fattibile
	double nsteps_tmp = t / stepsize;
	if(nsteps_tmp > (double) INT_MAX){
		fprintf(stderr, "Error: intregration step is too small\n");
		return -1;
	}
	
	int nsteps = (int) nsteps_tmp;
	
	// Alloco il vettore di output in base alla dimensione
	// del problema
	*output = malloc(nsteps * sizeof(double));
	if(*output == NULL) {
		fprintf(stderr, "Error: system memory is unavailable\n");
		return -1;
	}
	
	// Metodo di Eulero:
	(*output)[0] = x0;
	for(int i = 1; i < nsteps; i++) {
		double step = (*fnc_ptr)(i * stepsize, (*output)[i - 1]);
		(*output)[i] = (*output)[i - 1] + stepsize * step; // copio l'ouput nel vettore di output
	}
	
	// ritorno il numero di step del problema
	return nsteps;
}
