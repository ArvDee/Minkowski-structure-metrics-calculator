#include "main.h"

/* --------- ======= MAIN ======= ---------
 *
 */
int main(int argc, char *argv[]){
  if(argc < 3){
    printf("Incorrect number of arguments: please try './boop [MAX L] [FILE 1] [FILE 2] ... [FILE N]'\n");
    exit(42);
  }
  size_t max_l = strtol(argv[1],NULL,10);

  // Create a processor instance that will handle the calculation
  SnapshotProcessor processor;

  // Process all input files
  for(int arg = 2; arg < argc; arg++){
    // Check whether we have a coordinate text file (.dat, usually) or a GSD binary file (.gsd)
    std::string file_name(argv[arg]);
    // Strip the path from the file name
    std::string::size_type slash_idx = file_name.rfind('/');
    std::string short_file_name = file_name.substr(slash_idx+1, file_name.length());
    // Find the last dot in the file name
    std::string extension;
    std::string::size_type dot_idx = short_file_name.rfind('.');
    if(dot_idx != std::string::npos){
      extension = short_file_name.substr(dot_idx + 1);
    }else{
      printf("Error: No file extension found for supposed file '%s'.\n",file_name.c_str());
    }
    // Strip the extension from the file name and use short name to name BOOP files later
    short_file_name = short_file_name.substr(0, dot_idx);

    // Load the file and calculate bond order parameters
    if(extension == "gsd" || extension == "GSD"){
      // Create a GSD loader object that will handle the loading
      GSD_Loader gsd_loader;
      // Give the gsd loader pointers to where it should store the loaded data
    	gsd_loader.set_data_pointers(&processor.box, &processor.positions);
      // Open the file
      gsd_loader.open_gsd_file(file_name.c_str());
      size_t n_frames = gsd_loader.n_frames();
      // GSD files can contain multiple frames, so those should be specified by the user as well
      // The character '-' stands for "do all frames", "-1" for the last frame.
      std::string frame_str = "-1";
      if(arg+1 < argc) frame_str = std::string(argv[arg+1]);
      if( !strcmp(frame_str.c_str(),"-") ){ // strcmp returns 0 if string ARE identical, so invert
        for(size_t i = 0; i < n_frames; i++){
          gsd_loader.gsd_load_frame(i);
        }
      }else{
        char* non_nr;
        int frame_nr = strtol(frame_str.c_str(), &non_nr, 10);
        // Check to make sure the frame specifier is actually a number
        if(*non_nr && non_nr[0] != '-'){
          printf("Error: frame specifier '%s' is not '-' or a number.\n",frame_str.c_str());
          exit(42);
        }
        // Negative integers load backwards i.e. -1 is the last frame in the file.
        if(frame_nr < 0){
          frame_nr = n_frames - (frame_nr + 2); // +2 so that -1 becomes last element
        }
        // Load the specified frame
        gsd_loader.gsd_load_frame(frame_nr);
      }
      // The specified frame also takes up an argument slot => increment before processing next one
      arg++;
    }else if(extension == "dat"){
      // Load in the snapshot
      processor.load_snapshot(argv[arg]);
    }else{
      printf("Error: Could not recognize file extension '%s', try '.dat' or '.gsd'.\n",extension.c_str());
    }

    // Do the calculation
    processor.calculate_order_parameters(max_l);

    // Save the calculated bond order parameters to a file in the same directory as the source
    std::string dir_name = file_name.substr(0, slash_idx+1);
    processor.save_qw_files(dir_name,""); // second argument is an optional string to give the file name
  }
  
	return 1;
}
