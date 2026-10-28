add_library(physical_preparation
    physical_state.f90
    physical_particles.f90
    physical_preparation.f90
)
target_link_libraries(physical_preparation PUBLIC core)
