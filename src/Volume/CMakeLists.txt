file(GLOB SOURCES ./*.cpp )
file(GLOB HEADERS ./*.h)

add_library(volume_lib ${SOURCES} ${HEADERS})
target_link_libraries(volume_lib OpenVolumeMesh pmp_vis)

