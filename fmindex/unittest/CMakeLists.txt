include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/e2239ee6043f73722e7aa812a459f54a28552929.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Now simply link against gtest or gtest_main as needed. Eg
include_directories(${PROJECT_SOURCE_DIR}/src)
add_executable(week1 week1test.cpp ../src/fmindex.cpp )
add_executable(week2 week2test.cpp ../src/fmindex.cpp )
add_executable(week3 week3test.cpp ../src/fmindex.cpp ../src/bidirectionalfmindex.cpp )
target_link_libraries(week1 gtest_main)
target_link_libraries(week2 gtest_main)
target_link_libraries(week3 gtest_main)
