# ricpad library 
set(HEADER_PATH "include/ricpad")
set(HEADER 
    "${HEADER_PATH}/conversions.hpp"
    "${HEADER_PATH}/hankdet.hpp"
    )

add_library(ricpad STATIC src/conversions.cpp)
target_include_directories(ricpad PUBLIC include)
set_target_properties(ricpad PROPERTIES PUBLIC_HEADER "${HEADER}")
