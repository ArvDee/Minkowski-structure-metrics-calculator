#include "main.h"




unsigned int N_particles = 1;
double Tetrahedron_rounding = 0.0;
std::string initial_config_file_name;

/* --------- Process all command line input ---------
 * Sorts command line input into their relevant simulation components. Most notably, splits
 * simulation input from visualization input.
 */
void process_command_line_input(int argc, char *argv[]){
	 if(argc < 5){
 		printf("Incorrect number of arguments: please try ./fbmc -N [#particles] -s [ROUNDING]\n");
		printf("Optionally, add '-init [FILENAME]' to initialize from a certain configuration.\n");
 		exit(42);
 	}
 	if(argc >= 5){
 		for(int arg = 1; arg < argc; arg++){
 			if(strcmp(argv[arg],"-N") == 0){
 				N_particles = strtol(argv[arg+1],NULL,10);
 				if(N_particles == 0){printf("Invalid number of particles.\n");exit(42);}
 				printf("Argument '%s %s': set number of particles to %u.\n",argv[arg],argv[arg+1],N_particles);
 				arg+=1;
 			}
 			else if(strcmp(argv[arg],"-s") == 0){
 				Tetrahedron_rounding     = strtod(argv[arg+1],NULL);
 				printf("Argument '%s %s': set tetrahedra rounding fraction to %lf.\n",argv[arg],argv[arg+1],Tetrahedron_rounding);
 				if(Tetrahedron_rounding == 1.0){printf("Rounding is 1, making spheres.\n");} // Spheres
 				if(Tetrahedron_rounding == 0.0){printf("Rounding is 0, making sharp tetrahedra.\n");} // Sharp tetrahedra
				if(Tetrahedron_rounding < 0.0 || Tetrahedron_rounding > 1.0){
					printf("Please set a rounding between 0 and 1.\n");
					exit(42);
				}
 				arg+=1;
 			}
			else if(strcmp(argv[arg],"-init") == 0){
 				initial_config_file_name = std::string(argv[arg+1]);
 				printf("Argument '%s %s': loading initial configuration from file '%s'.\n",argv[arg],argv[arg+1],argv[arg+1]);
 				arg+=1;
 			}
			else{
				printf("Unknown input: '%s'. Please try './mc -N [#particles] -s [ROUNDING]'.\n",argv[arg]);
				exit(42);
			}
 		}
 	}
 }

/* --------- ======= MAIN ======= ---------
 *
 */
int main(int argc, char *argv[]){

	// Process command line input to split simulation and visualization input
	process_command_line_input(argc,argv);

	// Begin threaded environment where one thread does the simulation and possibly another the visualization

		// Perform simulation
		Simulation sim(N_particles, Tetrahedron_rounding, initial_config_file_name);

		// Also do live visualization
		// somethingsomething scene();

	// End threaded environment

	return 1; // Simulation ended successfully
}
