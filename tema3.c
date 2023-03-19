#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define min(a,b) (((a) < (b)) ? (a) : (b))

#define N_COORD 4 // numar procese coordonator

/* Numar si vector de vecini: pentru un proces coordonator vecinii vor fi 
workerii din cluster, iar pentru un proces de tip worker va exista un singur 
vecin, si anume coordonatorul */
static int num_neigh;
static int* neigh;

/* Verifica daca procesul dat este un coordonator: 1 - da, 0 - nu */
int is_coord(int rank) {
	return (0 <= rank) && (rank < N_COORD);
}

/* Pentru un coordonator citeste din fisierul corespunzator workerii din cluster, 
actualizand in vectorul de vecini; in cazul unui worker, aceasta se initializeaza 
cu dimensiune 1 */
void get_neighbours(int rank) {

	if (is_coord(rank)) {

		FILE *fp;
		char file_name[15];
		sprintf(file_name, "cluster%d.txt", rank);

		fp = fopen(file_name, "r");
		fscanf(fp, "%d", &num_neigh);

		neigh = malloc(sizeof(int) * num_neigh);

		for (size_t i = 0; i < num_neigh; i++)
			fscanf(fp, "%d", &neigh[i]);

	} else {
		num_neigh = 1;
		neigh = malloc(sizeof(int) * num_neigh);
	}
}

/* Logheaza un mesaj transmis de la sursa la destinatie */
void log_message(int src, int dest) {
	printf("M(%d,%d)\n", src, dest);
}

/* Pentru un coordonator dat, printeaza lista de workeri din cluster */
void print_workers(int coord, int* topology, int n_processes) {

	int first_print = 0;

	for (int i = 0; i < n_processes; i++) {
		if (topology[i] == coord) {
			if (first_print == 0) {
				printf("%d", i);
				first_print = 1;
			} else {
				printf(",%d", i);
			}
		}
	}
}

/* Afiseaza intreaga topologie, afisand lista de workeri pentru fiecare coordonator */
void print_topology(int* topology, int n_processes) {

	for (int i = 0; i < N_COORD; i++) {
		printf("%d:", i);
		print_workers(i, topology, n_processes);
		printf(" ");
	}
	printf("\n");
}

void gen_topology(int rank, int n_processes) {

	MPI_Status status;

	/* Topologie = vector de parinti */
	int* topology = malloc(sizeof(int) * n_processes);
	int* topology_recv = malloc(sizeof(int) * n_processes);

	memset(topology, -1, sizeof(int) * n_processes);
	memset(topology_recv, -1, sizeof(int) * n_processes);
	
	if (is_coord(rank)) {
		/* Un coordonator nu are parinte */
		topology[rank] = -1;

		/* Fiecare proces coordonator trimite un mesaj catre workerii sai pentru 
		ca acestia sa isi afle coordonatorul */
		for (int i = 0; i < num_neigh; i++) {
			MPI_Send(topology, n_processes, MPI_INT, neigh[i], 0, MPI_COMM_WORLD);
			log_message(rank, neigh[i]);
		}
	} else {
		/* Daca procesul curent este un worker, atunci va astepta un mesaj de 
		la coordonatorul sau, actualizand valoarea acestuia in vectorul de vecini 
		si in topologie */
		MPI_Recv(topology, n_processes, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		neigh[0] = topology[rank] = status.MPI_SOURCE;
	}

	if (!is_coord(rank)) {
		/* Un worker trimite topologia actualizata catre coordonatorul sau */
		MPI_Send(topology, n_processes, MPI_INT, neigh[0], 0, MPI_COMM_WORLD);
		log_message(rank, neigh[0]);

	} else {
		/* Procesul coordonator asteapta raspunsuri de la workeri, actualizand topologia */
		for (int i = 0; i < num_neigh; i++) {
			MPI_Recv(topology_recv, n_processes, MPI_INT, neigh[i], 0, MPI_COMM_WORLD, &status);

			for (int k = 0; k < n_processes; k++) {
				if (topology_recv[k] != -1) {
					topology[k] = topology_recv[k];
				}
			}
		}
	}

	/* Se trimite topologia de la 1 la 0 pentru a fi completata de toti coordonatorii
	(sens 1 -> 2 -> 3 -> 0) */

	if (is_coord(rank)) {

		/* Se calculeaza urmatorul rang si rangul precedent pe sensul 1 -> 2 -> 3 -> 0 */
		int next_rank = (rank == 3) ? 0 : (rank + 1);
		int prev_rank = (rank == 0) ? 3 : (rank - 1);

		if (rank == 1) {
			/* Procesul 1 trimite topologia catre 2 */
			MPI_Send(topology, n_processes, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
			log_message(rank, next_rank);

		} else {
			/* Se primeste topologia de la coordonatorul precedent, se actualizeaza topologia 
			si se trimite catre urmatorul coordonator */
			MPI_Recv(topology_recv, n_processes, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, &status);

			for (int k = 0; k < n_processes; k++) {
				if (topology_recv[k] != -1) {
					topology[k] = topology_recv[k];
				}
			}

			if (rank != 0) {
				/* Procesul 0, fiind ultimul, nu trebuie sa trimita mai departe topologia */
				MPI_Send(topology, n_processes, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
				log_message(rank, next_rank);
			}
		}
	}

	/* Topologia completa a ajuns la procesul 0, aceasta trebuie trimisa si celorlalti
	coordonatori care apoi o vor trimite catre workerii proprii (sens 0 -> 3 -> 2 -> 1) */

	if (is_coord(rank)) {

		/* Se calculeaza urmatorul rang si rangul precedent pe sensul 0 -> 3 -> 2 -> 1 */
		int next_rank = (rank == 0) ? 3 : (rank - 1);
		int prev_rank = (rank == 3) ? 0 : (rank + 1);

		if (rank == 0) {
			/* Procesul 0 trimite topologia catre 3 */
			MPI_Send(topology, n_processes, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
			log_message(rank, next_rank);

		} else {
			/* Se primeste topologia de la coordonatorul precedent si se trimite catre 
			urmatorul coordonator */
			MPI_Recv(topology, n_processes, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, &status);

			if (rank != 1) {
				/* Procesul 1, fiind ultimul, nu trebuie sa trimita mai departe topologia */
				MPI_Send(topology, n_processes, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
				log_message(rank, next_rank);
			}
		}

	}

	if (is_coord(rank)) {
		/* Fiecare coordonator trimite topologia completa catre workerii sai */
		for (int i = 0; i < num_neigh; i++) {
			MPI_Send(topology, n_processes, MPI_INT, neigh[i], 0, MPI_COMM_WORLD);
			log_message(rank, neigh[i]);
		}
	} else {
		/* Un worker primeste topologia completa de la coordonatorul sau */
		MPI_Recv(topology, n_processes, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
	}

	/* Se afiseaza topologia */
	printf("%d -> ", rank);
	print_topology(topology, n_processes);
}

void calculate_vector(int rank, int v_size, int n_processes) {

	MPI_Status status;

	int* v;
	int* v_recv;
	int* v_init;

	/* Initial se genereaza vectorul de procesul 0 si se trimite alaturi de dimensiunea
	acestuia de la 0 la 1, la toti coordonatorii (sens 0 -> 3 -> 2 -> 1) */

	if (is_coord(rank)) {

		/* Se calculeaza urmatorul rang si rangul precedent pe sensul 0 -> 3 -> 2 -> 1 */
		int next_rank = (rank == 0) ? 3 : (rank - 1);
		int prev_rank = (rank == 3) ? 0 : (rank + 1);

		if (rank == 0) {

			/* Se genereaza vectorul initial */
			v = malloc(sizeof(int) * v_size);

			for (int k = 0; k < v_size; k++) {
				v[k] = v_size - k - 1;
			}

			/* Se trimite vectorul initial si dimensiunea sa catre procesul 3 */
			MPI_Send(&v_size, 1, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
			MPI_Send(v, v_size, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
			log_message(rank, next_rank);

		} else {

			/* Asteapta dimensiunea vectorului de la coordonatorul precedent */
			MPI_Recv(&v_size, 1, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, &status);

			v = malloc(sizeof(int) * v_size);

			/* Asteapta vectorul de la coordonatorul precedent */
			MPI_Recv(v, v_size, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, &status);

			if (rank != 1) {
				/* Procesul 1, fiind ultimul, nu trebuie sa trimita mai departe topologia */
				MPI_Send(&v_size, 1, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
				MPI_Send(v, v_size, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
				log_message(rank, next_rank);
			}

		}

	}

	if (is_coord(rank)) {

		v_recv = malloc(sizeof(int) * v_size);

		/* Retine vectorul initial pentru a aplica ulterior modificarile realizate de workeri */
		v_init = malloc(sizeof(int) * v_size);
		for (int i = 0; i < v_size; i++) {
			v_init[i] = v[i];
		}

		/* Trimite vectorul alaturi de dimensiunea sa catre workerii sai, si limitele 
		intervalului din vector unde se vor realiza operatii */

		for (int i = 0; i < num_neigh; i++) {
			MPI_Send(&v_size, 1, MPI_INT, neigh[i], 0, MPI_COMM_WORLD);
			MPI_Send(v, v_size, MPI_INT, neigh[i], 0, MPI_COMM_WORLD);

			/* Imparte calculele in mod echilibrat */
			int start = (neigh[i] - N_COORD) * (double)v_size / (n_processes - N_COORD);
			int end = min((neigh[i] - N_COORD + 1) * (double)v_size / (n_processes - N_COORD), 
						v_size);

			MPI_Send(&start, 1, MPI_INT, neigh[i], 0, MPI_COMM_WORLD);
			MPI_Send(&end, 1, MPI_INT, neigh[i], 0, MPI_COMM_WORLD);

			log_message(rank, neigh[i]);

			/* Asteapta vectorul actualizat de la workeri si actualizeaza in coordonator */
			MPI_Recv(v_recv, v_size, MPI_INT, neigh[i], 0, MPI_COMM_WORLD, &status);

			for (int i = 0; i < v_size; i++) {
				if (v_recv[i] != v[i]) {
					v[i] = v_recv[i];
				}
			}

		}

	} else {

		/* Un worker primeste vectorul si dimensiunea sa de la coordonatorul sau */

		MPI_Recv(&v_size, 1, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);

		v = malloc(sizeof(int) * v_size);

		MPI_Recv(v, v_size, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);

		/* Primeste si intervalele vectorului unde se vor realiza operatii */

		int start, end;
		MPI_Recv(&start, 1, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&end, 1, MPI_INT, neigh[0], 0, MPI_COMM_WORLD, &status);

		/* Aplica operatia asupra elementelor din intervalul dat */
		for (int i = start; i < end; i++) {
			v[i] = v[i] * 5;
		}

		/* Trimite vectorul actualizat la coordonator */
		MPI_Send(v, v_size, MPI_INT, neigh[0], 0, MPI_COMM_WORLD);
		log_message(rank, neigh[0]);
	}

	/* Vectorul se trimite de la 1 la 0 pentru asamblarea acestuia (sens 1 -> 2 -> 3 -> 0) */

	if (is_coord(rank)) {

		/* Se calculeaza urmatorul rang si rangul precedent pe sensul 1 -> 2 -> 3 -> 0 */
		int next_rank = (rank == 3) ? 0 : (rank + 1);
		int prev_rank = (rank == 0) ? 3 : (rank - 1);

		if (rank == 1) {
			/* Se trimite vectorul catre procesul 2 */
			MPI_Send(v, v_size, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
			log_message(rank, next_rank);

		} else {

			/* Se primeste vectorul de la coordonatorul precedent, se actualizeaza  
			si se trimite catre urmatorul coordonator */
			MPI_Recv(v_recv, v_size, MPI_INT, prev_rank, 0, MPI_COMM_WORLD, &status);

			for (int i = 0; i < v_size; i++) {
				if (v_recv[i] != v_init[i]) {
					v[i] = v_recv[i];
				}
			}

			if (rank == 0) {

				/* Procesul 0 va contine vectorul complet si il va afisa */
				printf("Rezultat: ");
				for (int i = 0; i < v_size; i++) {
					printf("%d ", v[i]);
				}
				printf("\n");

			} else {
				/* Se trimite vectorul actualizat catre urmatorul coordonator */
				MPI_Send(v, v_size, MPI_INT, next_rank, 0, MPI_COMM_WORLD);
				log_message(rank, next_rank);
			}

		}

	}

}

int main(int argc, char *argv[]) {

	int rank, n_processes, v_size;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &n_processes);

	if (rank == 0) {
		v_size = atoi(argv[1]); // procesul 0 citeste dimensiunea vectorului
	}
	
	get_neighbours(rank);
	
	MPI_Barrier(MPI_COMM_WORLD);

	gen_topology(rank, n_processes);

	MPI_Barrier(MPI_COMM_WORLD);

	calculate_vector(rank, v_size, n_processes);
	
	MPI_Finalize();

	return 0;
}
