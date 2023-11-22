#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "simulate.h"

// Constants
#define SPATIAL_IMPACT 0.15

// Global variables
int num_points, num_steps, num_threads_global;
double *prev_wave, *current_wave, *next_wave, *final_wave;

// Function to compute one time step for a range of indices
void simulate_wave() {
#pragma omp parallel num_threads(num_threads_global)
  {
    int thread_id = omp_get_thread_num();
    int chunk_size = num_points / num_threads_global;
    int start_index = thread_id * chunk_size;
    int end_index = (thread_id == num_threads_global - 1)
                        ? num_points
                        : (thread_id + 1) * chunk_size;

    for (int t = 1; t <= num_steps; ++t) {
// Update each spatial point independently
#pragma omp for
      for (int i = start_index; i < end_index; ++i) {
        // Handle boundary conditions
        if (i != 0 && i != num_points - 1) {
          next_wave[i] =
              2 * current_wave[i] - prev_wave[i] +
              SPATIAL_IMPACT *
                  (prev_wave[i - 1] - 2 * current_wave[i] + prev_wave[i + 1]);
        }
      }

// Synchronize threads using barrier
#pragma omp barrier

// Copy the final results into the separate array
#pragma omp for
      for (int i = start_index; i < end_index; ++i) {
        final_wave[i] = next_wave[i];
      }

// Synchronize threads using barrier
#pragma omp barrier
    }
  }
}

// Main simulation function
double *simulate(const int i_max, const int t_max, const int num_threads,
                 double *old_array, double *current_array, double *next_array) {
  num_points = i_max;
  num_steps = t_max;
  num_threads_global = num_threads;
  prev_wave = old_array;
  current_wave = current_array;
  next_wave = next_array;
  final_wave = (double *)malloc(num_points * sizeof(double));

  // Start timer
  struct timeval start, end;
  gettimeofday(&start, NULL);

  // Run the simulation
  simulate_wave();

  // Copy the final results into the current_wave array
  for (int i = 0; i < num_points; ++i) {
    current_wave[i] = final_wave[i];
  }

  // Stop timer
  gettimeofday(&end, NULL);
  double execution_time =
      (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1e6;

  // Print execution time and normalized time
  printf("Wave Simulation Time: %f seconds\n", execution_time);

  if ((num_points * num_steps) > 0) {
    double normalized_time = execution_time / (num_points * num_steps);
    printf("Normalized Time: %f seconds\n", normalized_time);
  } else {
    printf("Error: Division by zero in normalized time calculation.\n");
  }

  free(final_wave);
  return current_wave;
}
