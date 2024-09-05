cmake_minimum_required(VERSION 3.16)

project(cnn_ckks_project)
set(CMAKE_CXX_COMPILER g++)

option(COMP "Build comparison operation" OFF)
option(BOOT "Build comparison operation" OFF)
option(CNN "Build comparison operation" ON)
option(AES_CTR "Build comparison operation" OFF)
option(GIFT64_CTR "Build comparison operation" OFF)

option(CHI2 "Build comparison operation" OFF)

# option(TEST "Build comparison operation" OFF)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./run)
find_package(SEAL 3.6 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(HEXL 1.2.5 HINTS ${OpenFHE_INTEL_HEXL_HINT_DIR} REQUIRED)
set(CMAKE_CXX_COMPILER g++-11)

if(TEST)
	add_subdirectory(common)
	add_subdirectory(cpu-ckks/single-key/comp)
	add_executable(gen_approx)
	target_sources(gen_approx
		PRIVATE
			./run/gen_approx.cpp
			${COMMON_SOURCE_FILES}
			${COMP_SOURCE_FILES}
	)
	target_compile_options(gen_approx PRIVATE -Wall -O3 -std=c++17 -pthread -march=native)
	target_link_libraries(gen_approx PRIVATE ntl gmp m pthread SEAL::seal)
	target_include_directories(gen_approx PUBLIC /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${COMP_HEADER_DIR})
endif()

if(COMP)
	add_subdirectory(common)
	add_subdirectory(cpu-ckks/single-key/comp)
	add_executable(compare)
	target_sources(compare
		PRIVATE
			./run/run_compare.cpp
			${COMMON_SOURCE_FILES}
			${COMP_SOURCE_FILES}
	)
	target_compile_options(compare PRIVATE -Wall -O3 -std=c++17 -pthread -march=native)
	target_link_libraries(compare PRIVATE ntl gmp m pthread SEAL::seal)
	target_include_directories(compare PUBLIC /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${COMP_HEADER_DIR})
endif()

if(BOOT)
	add_subdirectory(common)
	add_subdirectory(cpu-ckks/single-key/comp)
	if(NOT TARGET bootstrapping)
		add_subdirectory(cpu-ckks/single-key/ckks_bootstrapping)
	endif()
	add_executable(bootstrapping)
	target_sources(bootstrapping
		PRIVATE
			./run/run_bootstrapping.cpp
			${COMMON_SOURCE_FILES}
			${BOOTSTRAPPING_SOURCE_FILES}
	)
	target_compile_options(bootstrapping PRIVATE -Wall -O3 -std=c++17 -pthread -march=native)
	target_link_libraries(bootstrapping PRIVATE ntl gmp m pthread SEAL::seal)
	target_include_directories(bootstrapping PUBLIC /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${BOOTSTRAPPING_HEADER_DIR})
endif()


get_target_property(INTEL_HEXL_LOCATION HEXL::hexl LOCATION)
message(STATUS "INTEL_HEXL_LOCATION: ${INTEL_HEXL_LOCATION}")

get_target_property(INTEL_HEXL_INCLUDE_DIR HEXL::hexl INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "INTEL_HEXL_INCLUDE_DIR: ${INTEL_HEXL_INCLUDE_DIR}")

include_directories(${INTEL_HEXL_INCLUDE_DIR})


if(AES_CTR)
	add_subdirectory(common)
	add_subdirectory(cpu-ckks/single-key/seal_cipher)

	if(NOT TARGET compare)
		add_subdirectory(cpu-ckks/single-key/comp)
	endif()
	if(NOT TARGET bootstrapping)
		add_subdirectory(cpu-ckks/single-key/ckks_bootstrapping)
	endif()
	if(NOT TARGET aes_enc)
		add_subdirectory(cpu-ckks/single-key/aes_enc)
	endif()
	
	add_executable(aes_ctr)
	target_sources(aes_ctr
		PRIVATE
			./run/run_aes_ctr.cpp
			${COMMON_SOURCE_FILES}
			${SEAL_CIPHER_SOURCE_FILES}
			${AES_ENC_SOURCE_FILES}
			# ${COMP_SOURCE_FILES}
			${BOOTSTRAPPING_SOURCE_FILES}
			${INTEL_HEXL_INCLUDE_DIR}
	)
	target_compile_options(aes_ctr PRIVATE -w -O3 -std=c++17 -march=native -fopenmp)
	# target_link_libraries(aes_ctr PRIVATE ntl gmp m SEAL::seal)

	target_link_libraries(aes_ctr PRIVATE ntl gmp m SEAL::seal HEXL::hexl)
	target_include_directories(aes_ctr PUBLIC /usr/include /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${SEAL_CIPHER_HEADER_DIR} ${AES_ENC_HEADER_DIR} ${COMP_HEADER_DIR} ${BOOTSTRAPPING_HEADER_DIR} ${INTEL_HEXL_INCLUDE_DIR})
	# target_include_directories(cnn PUBLIC /usr/include /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${CNN_HEADER_DIR} ${COMP_HEADER_DIR})
endif()



if(GIFT64_CTR)
	add_subdirectory(common)
	add_subdirectory(cpu-ckks/single-key/seal_cipher)

	if(NOT TARGET compare)
		add_subdirectory(cpu-ckks/single-key/comp)
	endif()
	if(NOT TARGET bootstrapping)
		add_subdirectory(cpu-ckks/single-key/ckks_bootstrapping)
	endif()
	if(NOT TARGET gift64_ctr)
		add_subdirectory(cpu-ckks/single-key/gift64_enc)
	endif()
	
	add_executable(gift64_ctr)
	target_sources(gift64_ctr
		PRIVATE
			./run/run_gift64_ctr.cpp
			${COMMON_SOURCE_FILES}
			${SEAL_CIPHER_SOURCE_FILES}
			${GIFT64_ENC_SOURCE_FILES}
			# ${COMP_SOURCE_FILES}
			${BOOTSTRAPPING_SOURCE_FILES}
			${INTEL_HEXL_INCLUDE_DIR}
	)
	target_compile_options(gift64_ctr PRIVATE -w -O3 -std=c++17 -march=native -fopenmp)
	# target_link_libraries(aes_ctr PRIVATE ntl gmp m SEAL::seal)

	target_link_libraries(gift64_ctr PRIVATE ntl gmp m SEAL::seal HEXL::hexl)
	target_include_directories(gift64_ctr PUBLIC /usr/include /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${SEAL_CIPHER_HEADER_DIR} ${GIFT64_ENC_HEADER_DIR} ${COMP_HEADER_DIR} ${BOOTSTRAPPING_HEADER_DIR} ${INTEL_HEXL_INCLUDE_DIR})
	# target_include_directories(cnn PUBLIC /usr/include /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${CNN_HEADER_DIR} ${COMP_HEADER_DIR})
endif()

if(CHI2)
	add_subdirectory(common)
	add_subdirectory(cpu-ckks/single-key/seal_cipher)

	if(NOT TARGET compare)
		add_subdirectory(cpu-ckks/single-key/comp)
	endif()
	if(NOT TARGET bootstrapping)
		add_subdirectory(cpu-ckks/single-key/ckks_bootstrapping)
	endif()
	if(NOT TARGET aes_enc)
		add_subdirectory(cpu-ckks/single-key/aes_enc)
	endif()
	
	add_executable(chi2)
	target_sources(chi2
		PRIVATE
			./run/demo-chi2.cpp
			${COMMON_SOURCE_FILES}
			${SEAL_CIPHER_SOURCE_FILES}
			${AES_ENC_SOURCE_FILES}
			# ${COMP_SOURCE_FILES}
			${BOOTSTRAPPING_SOURCE_FILES}
			${INTEL_HEXL_INCLUDE_DIR}
	)
	target_compile_options(chi2 PRIVATE -w -O3 -std=c++17 -march=native -fopenmp)
	# target_link_libraries(aes_ctr PRIVATE ntl gmp m SEAL::seal)

	target_link_libraries(chi2 PRIVATE ntl gmp m SEAL::seal HEXL::hexl)
	target_include_directories(chi2 PUBLIC /usr/include /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${SEAL_CIPHER_HEADER_DIR} ${AES_ENC_HEADER_DIR} ${COMP_HEADER_DIR} ${BOOTSTRAPPING_HEADER_DIR} ${INTEL_HEXL_INCLUDE_DIR})
	# target_include_directories(cnn PUBLIC /usr/include /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${CNN_HEADER_DIR} ${COMP_HEADER_DIR})
endif()


if(CNN)
	add_subdirectory(common)
	add_subdirectory(cpu-ckks/single-key/seal_cipher)
	if(NOT TARGET compare)
		add_subdirectory(cpu-ckks/single-key/comp)
	endif()
	if(NOT TARGET bootstrapping)
		add_subdirectory(cpu-ckks/single-key/ckks_bootstrapping)
	endif()
	if(NOT TARGET cnn)
		add_subdirectory(cpu-ckks/single-key/cnn)
	endif()
	if(NOT TARGET aes_enc)
		add_subdirectory(cpu-ckks/single-key/aes_enc)
	endif()
	
	add_executable(cnn)
	target_sources(cnn
		PRIVATE
			./run/run_cnn.cpp
			${COMMON_SOURCE_FILES}
			${CNN_SOURCE_FILES}
			${SEAL_CIPHER_SOURCE_FILES}
			${AES_ENC_SOURCE_FILES}
			${COMP_SOURCE_FILES}
			${BOOTSTRAPPING_SOURCE_FILES}

	)
	target_compile_options(cnn PRIVATE -w -O3 -std=c++17 -pthread -march=native -fopenmp)
	target_link_libraries(cnn PRIVATE ntl gmp m pthread SEAL::seal gomp)
	target_include_directories(cnn PUBLIC /usr/include /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${SEAL_CIPHER_HEADER_DIR} ${AES_ENC_HEADER_DIR} ${CNN_HEADER_DIR} ${COMP_HEADER_DIR} ${BOOTSTRAPPING_HEADER_DIR})
	# target_include_directories(cnn PUBLIC /usr/include /usr/local/include /usr/local/include/NTL ${COMMON_HEADER_DIR} ${CNN_HEADER_DIR} ${COMP_HEADER_DIR})
endif()