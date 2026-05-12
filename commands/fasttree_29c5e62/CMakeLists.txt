cmake_minimum_required(VERSION 3.10)
project(FastTree C)

# Set the C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# --- Options ---
option(USE_SINGLE "Use single precision floating point (default is double)" OFF)
option(NO_SSE "Turn off use of SSE3 instructions" OFF)
option(USE_NATIVE "Optimize for the local machine architecture (-march=native)" ON)
option(USE_UNSAFE_MATH "Use -funsafe-math-optimizations" ON)
option(TRACK_MEMORY "Enable detailed reports of memory usage (slower)" OFF)

# --- Dependencies ---
# Math library
find_library(MATH_LIBRARY m)

# OpenMP
find_package(OpenMP)

# --- Compiler Flags ---
# Recommended flags from documentation
set(BASE_FLAGS -O0 -Wall)

if(USE_NATIVE)
    list(APPEND BASE_FLAGS -march=native)
endif()

if(USE_UNSAFE_MATH)
    list(APPEND BASE_FLAGS -funsafe-math-optimizations)
endif()

# Check for OpenMP SIMD support
set(SIMD_FLAGS -fopenmp-simd)

# Definitions based on options
if(USE_SINGLE)
    add_definitions(-DUSE_SINGLE)
endif()

if(NO_SSE)
    add_definitions(-DNO_SSE)
endif()

if(TRACK_MEMORY)
    add_definitions(-DTRACK_MEMORY)
endif()

# --- Source Files ---
set(FASTTREE_SOURCES 
        src/FastTree.c 
        src/FastNJ.c 
        src/hyper_parameters.c 
        src/command_line.c
        src/ReadAATransitionMatrix.c
        src/simd_vectors.c
        src/rev_functions.c
        src/trans_matrix.c
        src/random_number_generator.c
        src/TopHitJoin.c
        src/UniqueBestHits.c
        src/TopHitNJSearch.c
        src/hit_supports.c
        src/GetBestFromTopHits.c
        src/SetAllLeafTopHits.c
        src/abs_accuracy.c
        src/SetOutDistance.c
        src/CreateGTR.c
        src/AllocRateCategories.c
        src/ExhaustiveNJSearch.c
        src/FastNJSearch.c
        src/Free.c
        src/FreeNJ.c
        src/InitNJ.c
        src/Log.c
        src/Print.c
        src/ProgressReport.c
        src/ReadAlignment.c
        src/ReadTreeToken.c
        src/read_tree.c
        src/AddToFreq.c
        src/AlnToConstraints.c
        src/CorrectedPairDistances.c
        src/NormalizeFreq.c
        src/OutProfile.c
        src/ProfileDist.c
        src/ProfileDistPiece.c
        src/SeqDist.c
        src/SeqToProfile.c
        src/SetProfile.c
        src/UpdateOutProfile.c
        src/constraint.c
        src/ChooseNNI.c
        src/ExpEigenRates.c
        src/GTRNegLogLk.c
        src/MLQuartetLogLk.c
        src/MLQuartetNNI.c
        src/MLQuartetOptimize.c
        src/MLSiteLikelihoodsByRate.c
        src/MLSiteRates.c
        src/PDiffVector.c
        src/PSameVector.c
        src/PairLogLk.c
        src/PairNegLogLk.c
        src/PosteriorProfile.c
        src/ReadDistanceMatrix.c
        src/ReadMatrix.c
        src/ReadVector.c
        src/SetBestHit.c
        src/SetCodeDist.c
        src/SetMLGtr.c
        src/SetMLRates.c
        src/SetupDistanceMatrix.c
        src/TreeLength.c
        src/TreeLogLk.c
        src/ChooseNNI.c
        src/ExpEigenRates.c
        src/GTRNegLogLk.c
        src/MLQuartetLogLk.c
        src/MLQuartetNNI.c
        src/MLQuartetOptimize.c
        src/MLSiteLikelihoodsByRate.c
        src/MLSiteRates.c
        src/PDiffVector.c
        src/PSameVector.c
        src/PairLogLk.c
        src/PairNegLogLk.c
        src/PosteriorProfile.c
        src/ReadDistanceMatrix.c
        src/ReadMatrix.c
        src/ReadVector.c
        src/SetBestHit.c
        src/SetCodeDist.c
        src/SetMLGtr.c
        src/SetMLRates.c
        src/SetupDistanceMatrix.c
        src/TreeLength.c
        src/TreeLogLk.c
        src/OpenMPString.c
        src/MakeHashtable.c
        src/FindMatch.c
        src/hash.c
        src/clockDiff.c
        src/InitNNIStats.c
        src/NNI.c
        src/SPR.c
        src/UniquifyAln.c
        src/FreeNNIStats.c
        src/UpdateBranchLengths.c
        src/RecomputeProfiles.c
        src/OptimizeAllBranchLengths.c
        src/TestSplitsML.c
        src/RescaleGammaLogLk.c
        src/TestSplitsMinEvo.c
        src/ReliabilityNJ.c
        src/FreeUniquify.c
        src/FreeHashtable.c
        src/ResetTopVisible.c
        src/SetDistCriterion.c
        src/FreeProfile.c
        src/SetCriterion.c
        src/ActiveAncestor.c
        src/UpdateTopVisible.c
        src/UpdateVisible.c
        src/GetVisible.c
        src/NGaps.c
        src/gamma.c
        src/InitTraversal.c
        src/TraversePostorder.c
        src/FreeTraversal.c
        src/NewProfile.c
        src/brent_method.c
        src/RecomputeMLProfiles.c
        src/SetupABCD.c
        src/SkipTraversalInto.c
        src/DeleteUpProfile.c
        src/RecomputeProfile.c
        src/ReplaceChild.c
        src/UpdateForNNI.c
        src/FreeUpProfiles.c
        src/Sibling.c
        src/RootSiblings.c
        src/FindSPRSteps.c
        src/UnwindSPRStep.c
        src/UpProfiles.c
        src/GetUpProfile.c
        src/MLPairOptimize.c
        src/GammaLogLk.c
        src/ResampleColumns.c
        src/SHSupport.c
        src/SplitSupport.c
        src/QuartetWeight.c
        src/PathToRoot.c
        src/QuartetWeight.c
        src/FreePath.c
        src/main.c
        src/ReportSetting.c
        src/MakeMatrix.c
        src/ReadConstraints.c
        src/PrintStats.c
)

# --- Targets ---

# 1. FastTree (Serial)
add_executable(FastTree ${FASTTREE_SOURCES})
target_compile_options(FastTree PRIVATE ${BASE_FLAGS} ${SIMD_FLAGS})
target_include_directories(FastTree PRIVATE src)
if(MATH_LIBRARY)
    target_link_libraries(FastTree PRIVATE ${MATH_LIBRARY})
endif()

# 2. FastTreeMP (Multi-threaded)
if(OpenMP_C_FOUND)
    add_executable(FastTreeMP ${FASTTREE_SOURCES})
    # FastTreeMP requires -DOPENMP and -fopenmp
    target_compile_options(FastTreeMP PRIVATE ${BASE_FLAGS} ${SIMD_FLAGS} -DOPENMP)
    target_include_directories(FastTreeMP PRIVATE src)
    target_link_libraries(FastTreeMP PRIVATE OpenMP::OpenMP_C)
    if(MATH_LIBRARY)
        target_link_libraries(FastTreeMP PRIVATE ${MATH_LIBRARY})
    endif()
    message(STATUS "FastTreeMP will be built with OpenMP support")
else()
    message(STATUS "OpenMP not found. FastTreeMP will not be built.")
endif()

# --- Installation ---
install(TARGETS FastTree DESTINATION bin)
if(TARGET FastTreeMP)
    install(TARGETS FastTreeMP DESTINATION bin)
endif()
