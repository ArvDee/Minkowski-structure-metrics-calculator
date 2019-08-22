#ifndef GSD_LOADER_H
#define GSD_LOADER_H

#include "gsd.h"

// TEMP
#include <iostream>


/* ============================== Description =================================
 * More info on GSD format, such as chunk names and data types:
 * (https://gsd.readthedocs.io/en/latest/index.html)
 *===========================================================================*/


class GSD_Loader {
private:
  // Variables to load data into
  std::vector<float> *a1_; // lattice vectors
  std::vector<float> *a2_;
  std::vector<float> *a3_;
  std::vector<std::vector<float>> *positions_;
  // GSD file variables
  gsd_handle handle_;
	size_t n_frames_;
  // GSD data chunk variables
  const gsd_index_entry* chunk_;
  size_t chunk_size_;
  char* raw_data_;
public:
  GSD_Loader();
  virtual ~GSD_Loader();
  size_t n_frames(void){ return n_frames_; }
  // Function to set the arrays to load data into
  void set_data_pointers(
    std::vector<float> *a1,
    std::vector<float> *a2,
    std::vector<float> *a3,
    std::vector<std::vector<float>> *positions
  );
  // Opens a GSD file, sets handle_, and n_frames_
  void open_gsd_file(const char *gsd_file_name);
  // Attempt to load a data chunk with a specific name, sets chunk_, chunk_size_ and raw_data
  const gsd_index_entry* gsd_load_chunk(uint64_t frame, const char* chunk_name);
  // Main bridge function that loads in GSD data into partViewer variables
  void gsd_load_frame(uint64_t frame);
};

inline GSD_Loader::GSD_Loader():
  n_frames_(0), chunk_(nullptr), chunk_size_(0), raw_data_(nullptr)
{}

inline GSD_Loader::~GSD_Loader(void){
  // Close gsd file after loading all relevant data and free memory
  gsd_close(&handle_);
  free(raw_data_);
}

// Sets the pointers to the data arrays where the box and position data should go
inline void GSD_Loader::set_data_pointers(
  std::vector<float> *a1,
  std::vector<float> *a2,
  std::vector<float> *a3,
  std::vector<std::vector<float>> *positions
){
  a1_ = a1;
  a2_ = a2;
  a3_ = a3;
  positions_ = positions;
}

// Opens a GSD file, sets handle_, and n_frames_
inline void GSD_Loader::open_gsd_file(const char *gsd_file_name){
  // Open the GSD file and fill the handle
  switch (gsd_open(&handle_, gsd_file_name, GSD_OPEN_READONLY)){
		case 0: break; // successfully opened file
		case -1: printf("Error opening GSD file: IO error.\n");                  break;
		case -2: printf("Error opening GSD file: file is not a GSD file.\n");    break;
		case -3: printf("Error opening GSD file: invalid GSD version.\n");       break;
		case -4: printf("Error opening GSD file: corrupt file.\n");              break;
		case -5: printf("Error opening GSD file: failed to allocate memory.\n"); break;
		default: printf("Error: unknown return value from gsd_open.\n");         break;
	}
  // Get the number of frames in this gsd file
	n_frames_ = gsd_get_nframes(&handle_);
  if(n_frames_ == 0){ n_frames_ = 1; }
}

// Attempt to load a data chunk with a specific name, sets chunk_, chunk_size_ and raw_data
inline const gsd_index_entry* GSD_Loader::gsd_load_chunk(uint64_t frame, const char* chunk_name){
  // First get a pointer to the relevant chunk (here, the number of particles in the frame).
  chunk_ = gsd_find_chunk(&handle_, frame, chunk_name);

  // Check whether the chunk was opened successfully (errors currently handled in gsd_load_frame)
  // if(chunk_ == NULL){ printf("Warning: could not find '%s' chunk.\n", chunk_name); return nullptr; }
  // NOTE: this check _has_ to be for NULL, not nullptr, because gsd does it that way!
  if(chunk_ == NULL){ return chunk_; }

  // Allocate an array to store the raw data before conversion to the right type
  chunk_size_ = chunk_->N * chunk_->M * gsd_sizeof_type((gsd_type) chunk_->type); // in bytes
  raw_data_ = (char*) realloc(raw_data_, chunk_size_);

  // Read raw data
  switch (gsd_read_chunk(&handle_, raw_data_, chunk_)){
    case 0: break; // successfully loaded data chunk
    case -1: printf("Error loading GSD chunk: IO error.\n");               break;
    case -2: printf("Error opening GSD chunk: invalid input.\n");          break;
    case -3: printf("Error opening GSD chunk: invalid data file.\n");      break;
    default: printf("Error: unknown return value from gsd_read_chunk.\n"); break;
  }
  return chunk_;
}

// Main bridge function that loads in GSD data
inline void GSD_Loader::gsd_load_frame(uint64_t frame){
  // Check for some common errors
  if(frame > n_frames_){
    printf("Error: frame %lu/%lu is out of range!\n",frame, n_frames_-1);
  }else{
    printf("Loading frame %lu/%lu.\n",frame, n_frames_-1);
  }

  // Clear any data from previous loads
  positions_->clear();

  // Load number of particles
  if(gsd_load_chunk(frame, "particles/N") == NULL){
    printf("Fatal error: could not find 'particles/N' chunk.\n");
    exit(42);
  }
  uint32_t n_particles = (reinterpret_cast<uint32_t*>(raw_data_))[0];

	// Load the box size. Hoomd outputs box dimensions as defined here:
	// https://gsd.readthedocs.io/en/stable/schema-hoomd.html#chunk-configuration/box
	// https://hoomd-blue.readthedocs.io/en/stable/box.html
  if(gsd_load_chunk(frame, "configuration/box") == NULL){
    printf("Fatal error: could not find 'configuration/box' chunk.\n");
    exit(42);
  }
	float* bx = reinterpret_cast<float*>(raw_data_);
  float Lx=bx[0], Ly=bx[1], Lz=bx[2], xy=bx[3], xz=bx[4], yz=bx[5];
  *a1_ = {Lx,    0,     0};
  *a2_ = {xy*Ly, Ly,    0};
  *a3_ = {xz*Lz, yz*Lz, Lz};

  // Positions, offset by the center of the box due to different origin convention
  if( gsd_load_chunk(frame, "particles/position") == NULL){
    printf("Fatal error: could not find 'particles/position' chunk.\n");
    exit(42);
  }
	float* pos = reinterpret_cast<float*>(raw_data_);
  for(size_t i = 0; i < 3*n_particles; i+=3){
    std::vector<float> position = { pos[i], pos[i+1], pos[i+2] };
    positions_->push_back(position);
  }
  // for (size_t i = 0; i < positions_->size(); i++) {
  // for (size_t i = 0; i < 10; i++) {
  //   std::cout << (*positions_)[i][0] << ' ' << (*positions_)[i][1] << ' ' << (*positions_)[i][2] << 'r'<< '\n';
  // }
}

#endif // end header guard
