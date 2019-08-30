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
 * l'indirizzo di un vettore di elementi. 
 *
 * esempio: double * a = malloc(5 * sizeof(double))
 * crea un'area di memoria  di dimensione 5 |__|__|__|__|__|
 * il puntatore a punta alla prima cella:
 *    a --> |__|__|__|__|__|
 *
 * l'elemento *a è quindi il valore nella prima cella:
 *   |*a|__|__|__|__|
 *
 * oppure *(a+1) o a[1] puntano alla seconda cella
 *   |__|a[1]|__|__|__|
 *
 * l'operazione *a si chiama dereferenziazione e indica l'operazione
 * di prendere il valore di un indirizzo. Altre operazioni sono:
 * - incremento di puntatore: a + 1 è l'indirizzo della seconda cella
 *     |__|__|__|__|__|
 *          ^ 
 *          | 
 *        a + 1
 *
 * - indirizzo di una variabile: &a ritorna l'indirizzo alla cella di memoria
 *   che contiene una variabile. Sia a = 5, &a è l'indirizzo della cella di memoria
 *   che contiene 5
 *
 * Ora l'utilizzo principale dei puntatori come argomento di funzione è quello di
 * modificare il valore di una variabile anche fuori dalla funzione. Ad esempio
 *    void function(int * a) { *a = 5; }
 * permette anche al chiamante di vedere la modifica a = 5. Questo perché alla funzione
 * viene passato un indirizzo di memoria e NON un valore. Andando a scrivere dentro
 * all'indirizzo, quando si esce dalla funzione si trova il valore modificato. Quindi
 * sia che facciamo:
 *   int a = 3;
 *   function(&a);
 *   printf("%d\n", a); // print 5
 *
 * sia che facciamo:
 *   int * a;
 *   *a = 3;
 *   function(a);
 *   printf("%d\n"); // print 5
 *
 * otteniamo lo stesso identico risultato. In entrambi i casi la funzione ha scritto 
 * in un'area di memoria puntata da un indirizzo, quindi quell'area viene modificata
 * per tutto il programma
 *
 * Il doppio puntatore serve perché noi non possiamo
 * conoscere da principio la dimensione dell'output del problema (che invece ritorniamo
 * come valore di ritorno della funzione) ma riallochiamo la dimensione del vettore
 * all'interno della funzione stessa.
 * Per far si che possiamo correttamente riallocare il vettore e rendere visibile la
 * modifica al di fuori della funzione ODE, dobbiamo modificare non solo il contenuto
 * del vettore, ma anche il suo indirizzo (quando allochiamo con malloc, il sistema
 * operativo ci ritorna un indirizzo all'area di memoria allocata). Da questo la
 * necessità del doppio puntatore.
 *
 * in questo caso:
 * double ** a è l'indirizzo all'indirizzo dell'array
 *
 *   a --> *a --> NULL
 *
 * ora passiamo alla malloc *a, cioè il puntatore all'area di memoria che vogliamo allocare
 * in questo modo:
 *   *a = malloc(...)
 *
 *   a --> *a --> |__|__|__|__| ... |__|
 *
 * abbiamo scritto nell'INDIRIZZO puntato da a l'INDIRIZZO di un vettore appena allocato.
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
	// del problema e lo scrivo nell'indirizzo di output
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


int ODE1_3D(double x0[3], double t, double stepsize, double **output, void(*fnc_ptr)(double, double *, double *))
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
	
	/* the output vector is allocated. The real size of the
	 * output is not 3 * # steps but is 4 * # steps. This will
	 * allow the use of 32byte vectorization at the cost of a
	 * greater memory usage. The wasted memory is 
	 * nsteps * sizeof(double). 
	 *
	 * NOTE: the ALIGN_SIZE macro can be set to 3 if 32BYTE_VECTORIZATION
	 * macro is disabled
	 */
	*output = malloc(ALIGN_SIZE * nsteps * sizeof(double));
	if(*output == NULL) {
		fprintf(stderr, "Error: system memory is unavailable\n");
		return -1;
	}
	
	double * step_tmp = malloc(4 * sizeof(double));
	if(step_tmp == NULL) {
		fprintf(stderr, "Error: system memory is unavailable\n");
		return -1;
	}
	
	// Memory initialization. In the *output vector the first three
	// points are set as the x0 values
	memcpy(*output, x0, 3 * sizeof(double));
	// then the step vector is initialized with 0s. Is important to
	// initialize the whole vector (expecially if its size is 4 due
	// to alignment purposes, to avoid getting random numbers in the
	// output)
	memset(step_tmp, 0, ALIGN_SIZE * sizeof(double));

	// vectorization helper pointers are created here.
	double * v_in;
	double * v_out;
	for(int i = 1; i < nsteps; i++) {
		
		// This vector manipulation is used to prepare vectorization.
		// It seems that gcc does not like working with double pointers
		// when trying to vectorize
		v_in = (*output) + nextindex(i - 1);
		v_out = v_in + ALIGN_SIZE;
		
		(*fnc_ptr)(i * stepsize, v_in, step_tmp);
		
		/* the unroll procedure is used to prepare vectorization.
		 * the auto-vectorization is requested to the compiler, so
		 * no explicit instructions are needed. However in order to
		 * perform a successfull SIMD set compilation, the for loop
		 * must be unrolled and the operations must be simple.
		 *
		 * With AVX and SSE instruction sets, only 2 doubles can be
		 * operated at the same time, this means that the two for loop 
		 * will be completed in 4 operations instead of 6
		 *
		 * with AVX2 and SSE3 instruction sets, up to 4 doubles can be
		 * operated at the same time, this means that the two for loop
		 * will be completed in 2 operations instead of 6.
		 *
		 * Moreover this allow us to work with quaternions in 2 operation
		 * as a future developement.
		 *
		 * Refer to the make file to enable SIMD compilation of the code
		 */
		#pragma unroll(4)
		#pragma simd
		for(int j = 0; j < ALIGN_SIZE; j++) {
			step_tmp[j] = step_tmp[j] * stepsize;
		}
		
		/* the output vector will be formatted in the following way:
		 * |t0_x0|t0_x1|t0_x2|EMPTY|t1_x0|t1_x1|t1_x2|EMPTY| ... |tn_x0|tn_x1|tn_x2|EMPTY|
		 * the EMPTY element is used to perform 32byte vector alignment.
		 * this will be used to perform 4 operations at time
		 */
		
		#pragma unroll(4)
		#pragma simd
		for(int j = 0; j < ALIGN_SIZE; j++) {
			v_out[j] = v_in[j] + step_tmp[j];
		}
	}
	// The return value is the number of steps performed (that is size of output / ALIGN_SIZE)
	return nsteps;
}
