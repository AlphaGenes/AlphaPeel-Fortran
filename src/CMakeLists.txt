SET(AlphaMLPEXE AlphaPeel)
SET(AlphaMLP_src AlphaMLPInputModule.f90
AlphaMLP.f90
pedigreeGraphModule.f90
secantEstimator.f90
globalMLP.f90
graphs.f90)

if(GPLIB)

    ADD_LIBRARY(${GPLIB} STATIC ${AlphaMLP_src})
    TARGET_LINK_LIBRARIES(${GPLIB} ${AHLIB})
    TARGET_LINK_LIBRARIES(${GPLIB} ${MKL_BLAS_LIB})
# we won't build exe and lib at the same time
else()
    ADD_EXECUTABLE(${AlphaMLPEXE}  ${AlphaMLP_src}  ${SRCAlphaMLP}main.f90)

    TARGET_LINK_LIBRARIES(${AlphaMLPEXE} ${AHLIB})
    #    TARGET_LINK_LIBRARIES(${AlphaMLPEXE} ${HMMLIB})
    TARGET_LINK_LIBRARIES(${AlphaMLPEXE} ${MKL_BLAS_LIB})
endif()



