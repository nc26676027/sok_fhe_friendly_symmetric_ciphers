#include <iostream>
#include <unistd.h>   // 对于 fork(), getpid()
#include <sys/wait.h> // 对于 waitpid()
#include <omp.h>      // 对于 OpenMP
#include <vector>
#include <bitset>
#include <chrono>


// #include "seal/seal.h"
#include "examples.h"
using namespace std;

// 定义一个简单的累加函数
bool getrandbit(int times) {
  bool tmp = 0;
  // If state has not been initialized yet
  std::bitset<80> state;

  if (state.none()) {
    state.set();  // Initialize with all bits set
    // Throw the first 160 bits away
    for (unsigned i = 0; i < 160; ++i) {
      // Update the state
      tmp = state[0] ^ state[13] ^ state[23] ^ state[38] ^ state[51] ^ state[62];
      state >>= 1;
      state[79] = tmp;
    }
  }
  // choice records whether the first bit is 1 or 0.
  // The second bit is produced if the first bit is 1.
  bool choice = false;
  do {
    // Update the state
    tmp = state[0] ^ state[13] ^ state[23] ^ state[38] ^ state[51] ^ state[62];
    state >>= 1;
    state[79] = tmp;
    choice = tmp;
    tmp = state[0] ^ state[13] ^ state[23] ^ state[38] ^ state[51] ^ state[62];
    state >>= 1;
    state[79] = tmp;
  } while (!choice);
//   cout << "in get randbit and run "<<times<<endl;
  return tmp;
}

void single_time(int threads){
    cout << "this point, in single time!"<<endl;
    #pragma omp parallel for num_threads(threads)
    for(int i=0;i<256;i++){
        // evaluator.multiply_inplace( mult_vec_1[i], mult_vec_2[i] );
        // evaluator.relinearize_inplace(mult_vec_1[i], relin_keys);
        // evaluator.rescale_to_next_inplace(mult_vec_1[i]);
        for(int j=0;j<327680;j++){
            bool my = getrandbit(j);
        }
    }
}


void test_threads(){
    pid_t pid = fork();
    

    if (pid == -1) {
        // fork 失败
        perror("fork");
    } else if (pid == 0) {
        // 子进程
        std::cout << "Child Process ID: " << getpid() << std::endl;
        auto child_st = std::chrono::high_resolution_clock::now();
        single_time(64);
        auto child_ed = std::chrono::high_resolution_clock::now();
        auto child_duration = std::chrono::duration_cast<std::chrono::milliseconds>(child_ed - child_st);
        std::cout <<4<<"Threads, child time taken: " << child_duration.count() << " milliseconds" << std::endl;
 
    } else {
        // 父进程
        std::cout << "Parent Process ID: " << getpid()
                  << ", Child PID: " << pid << std::endl;
        auto parent_st = std::chrono::high_resolution_clock::now();
        single_time(64);
        auto parent_ed = std::chrono::high_resolution_clock::now();
        auto parent_duration = std::chrono::duration_cast<std::chrono::milliseconds>(parent_ed - parent_st);
        std::cout <<4<<"Threads, parent time taken: " << parent_duration.count() << " milliseconds" << std::endl;
  
        // 等待子进程结束
        waitpid(pid, NULL, 0);
    }
}



void try_multi_threads() {

    
    auto start = std::chrono::high_resolution_clock::now();

    test_threads();
    // 获取结束时间点
    auto end = std::chrono::high_resolution_clock::now();
    // 计算所经过的时间，以毫秒为单位
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // 输出运行时间
    std::cout <<4<<"Threads, Time taken: " << duration.count() << " milliseconds" << std::endl;

}