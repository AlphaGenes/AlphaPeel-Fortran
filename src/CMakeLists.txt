SET(AlphaMLPEXE AlphaMLP)
SET(AlphaMLP_src ${SRCAlphaMLP}AlphaMLPInputModule.f90
 ${SRCAlphaMLP}AlphaMLP.f90
${SRCAlphaMLP}pedigreeGraphModule.f90
 ${SRCAlphaMLP}secantEstimator.f90
 ${SRCAlphaMLP}globalMLP.f90
 ${SRCAlphaMLP}graphs.f90
            ../${SRCALPHAIMPUTE}AlphaImpute.f90
            ../${SRCALPHAIMPUTE}AlphaImputeInMod.f90
            ../${SRCALPHAIMPUTE}global.f90
            ../${SRCALPHAIMPUTE}HaplotypeBits.f90
            ../${SRCALPHAIMPUTE}hmm.f90
            ../${SRCALPHAIMPUTE}Output.f90
            ../${SRCALPHAIMPUTE}hmmvariables.f90
            ../${SRCALPHAIMPUTE}hmmHaplotype.f90      
            ../${SRCALPHAIMPUTE}informationModule.f90    
            ../${SRCALPHAIMPUTE}Imputation.f90
            ../${SRCALPHAIMPUTE}InputMod.f90
            ../${SRCALPHAIMPUTE}manageWindows.f90
            ../${SRCALPHAIMPUTE}par_zig_mod.f90
            ../${SRCALPHAIMPUTE}PhaseRounds.f90
            ../${SRCALPHAIMPUTE}random.f90
            ../${SRCALPHAIMPUTE}utils.f90)

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



