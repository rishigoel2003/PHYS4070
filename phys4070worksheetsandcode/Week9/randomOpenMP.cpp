#include <iostream>
#include <omp.h> //for omp_get_thread_num() etc.
#include <random>
#include <utility> // for std::pair (not needed on most compilers)

// This function is used in the last example (method 4)
auto my_generator() {
  // static: will only be constructed ONCE, will live for life of program
  // nb: thread_local implies static: i.e., thread_local = static thread_local

  // std::random_device{}() is one-line version of:
  // std::random_device rd; rd();

  // I add omp_get_thread_num() to this (the seed) to ensure each thread
  // gets a different seed (but thread_local variables should be initialised in
  // thread-safe way, so this is strictly not required [I think]).

  thread_local std::mt19937 generator_tl(std::random_device{}() +
                                         omp_get_thread_num());

  return generator_tl();
}

//============================================================================
//============================================================================
int main() {

  std::cout << "===========================================\n";
  std::cout << "Random number generation with OpenMP:\n";
  std::cout << "===========================================\n";

  // Number of random numbers to generate
  int n_loop = 8;

  //============================================================================
  // METHOD 1: BAD, RACE CONDITION
  //============================================================================
  std::cout << "Random numbers in parallel, with race condition:\n\n";
  std::cout
      << "The race condition is caused by each thread using the same random\n"
         "number generator - there is a chance two threads will try to do "
         "this\n"
         "at the same time, so they will get the _same_ random number (bad!)\n";


  // list to store calculated random numbers, and which thread
  std::vector<std::pair<int, int>> list1(n_loop);

  // Initiate a Mersenne Twister generator, seed it using random_device
  int seed = std::random_device{}();
  std::mt19937 generator_1(seed);

  // Parallel for loop
  #pragma omp parallel for
  for (int i = 0; i < n_loop; ++i) {
    // Get thread number, random number, store in the list of pairs
    int thread = omp_get_thread_num();
    int rn = generator_1();
    list1.at(i) = {thread, rn};
  }
  // End parallel region

  // Serial for loop over vector "list1" - print out random numbers from each thread
  for (const std::pair<int, int> &item : list1) {
    int thread = item.first;
    int rand = item.second;
    std::cout << thread << ": " << rand << "\n";
  }
  std::cout << "\n";

  //============================================================================
  // METHOD 2: CLOSER, BUT STILL BAD
  //============================================================================
  std::cout << "---------------------------------------------------------\n\n";

  std::cout
      << "Random numbers in parallel, with thread_local (incorrectly):\n\n";
  std::cout
      << "We can avoid this race condition by declaring our generator to\n"
         "be `thread_local` - this means the variable will be stored\n"
         "seperately on each thread that initialises it.\n"
         "The problem with the way we do it here, is that only one thread "
         "initialises the generator, so the other threads don't have access\n";

  // list to store calculated random numbers, and which thread
  std::vector<std::pair<int, int>> list2(n_loop);

  // Generate a "thread local" random number generator
  // Problem: this is outside the parallel region, so it will only initiate
  //          the generator on thread 0 (the master thread)
  thread_local std::mt19937 generator_tl1(std::random_device{}() +
                                          omp_get_thread_num());
  // nb: thread_local implies static: i.e., thread_local = static thread_local

  // Parallel for loop
  #pragma omp parallel for
  for (int i = 0; i < n_loop; ++i) {
    // Get thread number, random number, store in the list of pairs
    int thread = omp_get_thread_num();
    int rn = generator_tl1();
    list2.at(i) = {thread, rn};
  }
  // End parallel region

  // Serial for loop over vector "list2" - print out random numbers from each thread
  for (const std::pair<int, int> &item : list2) {
    int thread = item.first;
    int rand = item.second;
    std::cout << thread << ": " << rand << "\n";
  }
  std::cout << "\n";

  //============================================================================
  // METHOD 3: INITIATE RANDOM NUMBER GENERATOR INSIDE PARALLEL REGION
  //============================================================================
  std::cout << "---------------------------------------------------------\n\n";
  std::cout << "Random numbers in parallel, with thread_local, method 1:\n\n";
  std::cout << "To correct this, we have to make this call _inside_ the "
               "parellel region:\n";

  // list to store calculated random numbers, and which thread
  std::vector<std::pair<int, int>> list3(n_loop);

  // Start parallel region
  #pragma omp parallel
  {
    // Generate a "thread local" random number generator
    thread_local std::mt19937 generator_tl(std::random_device{}() +
                                           omp_get_thread_num());
    // nb: thread_local implies static: i.e., thread_local = static thread_local

    // Parallel for loop
    #pragma omp for
    for (int i = 0; i < n_loop; ++i) {
      // Get thread number, random number, store in the list of pairs
      int thread = omp_get_thread_num();
      int rn = generator_tl();
      list3.at(i) = {thread, rn};
    }
  }
  // End parallel region

  // Serial for loop over vector "list3" - print out random numbers from each thread
  for (const std::pair<int, int> &item : list3) {
    int thread = item.first;
    int rand = item.second;
    std::cout << thread << ": " << rand << "\n";
  }
  std::cout << "\n";

  //============================================================================
  // METHOD 4: USE A THREAD_LOCAL FUNCTION (DEFINED AT THE TOP OF THE SCRIPT)
  //============================================================================
  std::cout << "---------------------------------------------------------\n\n";

  std::cout << "Random numbers in parallel, with thread_local in function:\n\n";
  std::cout
      << "Making the call inside the parallel region can be annoying and\n"
         "limitting. Instead, we can wrap this in a new function that\n"
         "contains a thread_local variable. Thread_local variables are\n"
         "also static, meaning they will only be initialised ONCE, and\n"
         "then live for the life of the program. thread_local means it\n"
         "will be initialised once PER THREAD\n";

  // list to store calculated random numbers, and which thread
  std::vector<std::pair<int, int>> list4(n_loop);

  // Parallel for loop
  #pragma omp parallel for
  for (int i = 0; i < n_loop; ++i) {
    int thread = omp_get_thread_num();
    int rn = my_generator();
    list4.at(i) = {thread, rn};
  }
  // End parallel region

  // Serial for loop over vector "list3" - print out random numbers from each thread
  for (const std::pair<int, int> &item : list4) {
    int thread = item.first;
    int rand = item.second;
    std::cout << thread << ": " << rand << "\n";
  }
}