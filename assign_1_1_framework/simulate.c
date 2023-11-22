/*
 * simulate.c
 *
 * Implement your (parallel) simulation here!
 */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "simulate.h"

#define SPATIAL_IMPACT 0.15

int num_points, num_steps, num_threads_global;
double *prev_wave, *current_wave, *next_wave, *final_wave;
pthread_mutex_t mutex;
pthread_cond_t cond;
int barrier_count = 0;

// Function to compute one time step for a range of indices
void *simulate_wave(void *arg) {
  long thread_id = (long)arg;
  int chunk_size = num_points / num_threads_global;
  int start_index = thread_id * chunk_size;
  int end_index = (thread_id == num_threads_global - 1)
                      ? num_points
                      : (thread_id + 1) * chunk_size;

  for (int t = 1; t <= num_steps; t++) {
    for (int i = start_index; i < end_index; i++) {
      if (i == 0 || i == num_points - 1) {
        next_wave[i] = current_wave[i];
      } else {
        next_wave[i] =
            2 * current_wave[i] - prev_wave[i] +
            SPATIAL_IMPACT *
                (prev_wave[i - 1] - 2 * current_wave[i] + prev_wave[i + 1]);
      }
    }

    // Synchronize threads using mutex and conditional variable
    pthread_mutex_lock(&mutex);
    barrier_count++;
    if (barrier_count == num_threads_global) {
      barrier_count = 0;
      pthread_cond_broadcast(&cond);
    } else {
      while (pthread_cond_wait(&cond, &mutex) != 0)
        ;
    }
    pthread_mutex_unlock(&mutex);
  }

  // Copy the final results into array
  for (int i = start_index; i < end_index; ++i) {
    final_wave[i] = next_wave[i];
  }

  return NULL;
}

/*
 * Executes the entire simulation.
 *
 * Implement your code here.
 *
 * i_max: how many data points are on a single wave
 * t_max: how many iterations the simulation should run
 * num_threads: how many threads to use (excluding the main threads)
 * old_array: array of size i_max filled with data for t-1
 * current_array: array of size i_max filled with data for t
 * next_array: array of size i_max. You should fill this with t+1
 */
double *simulate(const int i_max, const int t_max, const int num_threads,
                 double *old_array, double *current_array, double *next_array) {
  num_points = i_max;
  num_steps = t_max;
  num_threads_global = num_threads;
  prev_wave = old_array;
  current_wave = current_array;
  next_wave = next_array;

  // Initialize mutex and conditional variable
  pthread_mutex_init(&mutex, NULL);
  pthread_cond_init(&cond, NULL);

  // Allocate memory for the final results array
  final_wave = (double *)malloc(num_points * sizeof(double));

  // Create threads
  pthread_t threads[num_threads];
  for (long i = 0; i < num_threads; i++) {
    printf("numthreads: %d\n", num_threads);
    pthread_create(&threads[i], NULL, simulate_wave, (void *)i);
  }

  struct timeval start, end;
  gettimeofday(&start, NULL);

  // Main thread waits for worker threads to finish
  for (int t = 0; t < num_threads; ++t) {
    pthread_join(threads[t], NULL);
  }

  // Copy the final results into the current_wave array
  for (int i = 0; i < num_points; ++i) {
    current_wave[i] = final_wave[i];
  }

  gettimeofday(&end, NULL);
  double execution_time =
      (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec);

  pthread_mutex_destroy(&mutex);
  pthread_cond_destroy(&cond);

  printf("Wave Simulation Time: %f seconds\n", execution_time);
  double normalized_time = execution_time / (num_points * num_steps);
  printf("Normalized Time: %f seconds\n", normalized_time);
  free(final_wave);
  return current_wave;
}
