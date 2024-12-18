# CMake generated Testfile for 
# Source directory: /home/f/Documents/SpecialTopicsInFiniteElementsBlue/tests
# Build directory: /home/f/Documents/SpecialTopicsInFiniteElementsBlue/cmake-build-debug/tests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(tests/template.debug "/snap/clion/308/bin/cmake/linux/x64/bin/cmake" "-DTRGT=tests.template.debug.test" "-DTEST=tests/template.debug" "-DEXPECT=PASSED" "-DBINARY_DIR=/home/f/Documents/SpecialTopicsInFiniteElementsBlue/cmake-build-debug" "-P" "/usr/share/deal.ii//scripts/run_test.cmake")
set_tests_properties(tests/template.debug PROPERTIES  LABEL "tests" PROCESSORS "1" TIMEOUT "600" WORKING_DIRECTORY "/home/f/Documents/SpecialTopicsInFiniteElementsBlue/cmake-build-debug/tests/template.debug" _BACKTRACE_TRIPLES "/usr/share/deal.ii/macros/macro_deal_ii_add_test.cmake;577;add_test;/usr/share/deal.ii/macros/macro_deal_ii_pickup_tests.cmake;343;deal_ii_add_test;/home/f/Documents/SpecialTopicsInFiniteElementsBlue/tests/CMakeLists.txt;2;DEAL_II_PICKUP_TESTS;/home/f/Documents/SpecialTopicsInFiniteElementsBlue/tests/CMakeLists.txt;0;")
