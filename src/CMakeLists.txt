set(HEADERS Surface/Curvature.h  Surface/diffgeo.h Surface/GeodesicsInHeat.h  Surface/Parameterization.h Surface/Smoothing.h  Surface/SpectralProcessing.h Surface/LaplaceConstruction.h Surface/Diamond_2D.h Surface/Poisson_System.h Volume/LaplaceConstruction_3D.h Volume/Diamond_3D.h Volume/Eigenmodes.h Volume/Franke_PoissonSystem_3D.h Volume/diffgeo_3D.h Volume/GeodesicsInHeat_3d.h)
set(SOURCES Surface/Curvature.cpp  Surface/diffgeo.cpp Surface/GeodesicsInHeat.cpp  Surface/Parameterization.cpp Surface/Smoothing.cpp  Surface/SpectralProcessing.cpp Surface/LaplaceConstruction.cpp Surface/Diamond_2D.cpp Surface/Poisson_System.cpp Volume/LaplaceConstruction_3D.cpp Volume/Diamond_3D.cpp Volume/Eigenmodes.cpp Volume/Franke_PoissonSystem_3D.cpp Volume/diffgeo_3D.cpp Volume/GeodesicsInHeat_3d.cpp)

add_subdirectory(VolumeMesh)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/VolumeMesh)

if (EMSCRIPTEN)

    add_library(poly STATIC ${HEADERS} ${SOURCES})
    target_link_libraries(poly pmp_vis VolumeMesh)

    add_executable(diamond_laplacian main.cpp SurfaceViewer.cpp SurfaceViewer.h VolumeViewer.cpp VolumeViewer.h)
    target_link_libraries(diamond_laplacian poly VolumeMesh)

    set_target_properties(diamond_laplacian PROPERTIES LINK_FLAGS "--shell-file ${PROJECT_SOURCE_DIR}/external/pmp-library/src/apps/data/shell.html")

else()

    add_library(poly STATIC ${HEADERS} ${SOURCES})
    target_link_libraries(poly pmp VolumeMesh)

    add_executable(diamond_laplacian main.cpp SurfaceViewer.cpp SurfaceViewer.h VolumeViewer.cpp VolumeViewer.h)
    target_link_libraries(diamond_laplacian poly pmp_vis VolumeMesh)


endif()
