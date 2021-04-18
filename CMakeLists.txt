cmake_minimum_required (VERSION 3.6)


project(Atmosphere CXX)

set(DILIGENT_ENGINE_PATH "D:/indie/DiligentEngine")


include(BuildUtils.cmake)


set(SOURCE
    src/AtmosphereSample.cpp
    src/Terrain/EarthHemisphere.cpp
    src/Terrain/ElevationDataSource.cpp
	
)

set(INCLUDE
    src/AtmosphereSample.hpp
    src/Terrain/DynamicQuadTreeNode.hpp
    src/Terrain/EarthHemisphere.hpp
    src/Terrain/ElevationDataSource.hpp
    src/Terrain/HierarchyArray.hpp
)

set(TERRAIN_SHADERS
    assets/shaders/terrain/GenerateNormalMapPS.fx
    assets/shaders/terrain/HemispherePS.fx
    assets/shaders/terrain/HemisphereVS.fx
    assets/shaders/terrain/HemisphereZOnlyVS.fx
    assets/shaders/terrain/ScreenSizeQuadVS.fx
    assets/shaders/terrain/TerrainShadersCommon.fxh
)

set(SHADERS
    assets/shaders/HostSharedTerrainStructs.fxh
	../DiligentEngine/DiligentFX/Shaders/PostProcess/EpipolarLightScattering/public/EpipolarLightScatteringStructures.fxh
    ../DiligentEngine/DiligentFX/Shaders/PostProcess/EpipolarLightScattering/public/EpipolarLightScatteringFunctions.fxh
    ../DiligentEngine/DiligentFX/Shaders/PostProcess/ToneMapping/public/ToneMappingStructures.fxh
    ../DiligentEngine/DiligentFX/Shaders/Common/public/BasicStructures.fxh
    ../DiligentEngine/DiligentFX/Shaders/Common/public/Shadows.fxh
)

set(TERRAIN_TEXTURES
    assets/Terrain/HeightMap.tif
    assets/Terrain/Mask.png
)

set(TERRAIN_TILES
    assets/Terrain/Tiles/cliff_DM.dds
    assets/Terrain/Tiles/cliff_NM.dds
    assets/Terrain/Tiles/grass_DM.dds
    assets/Terrain/Tiles/grass_NM.dds
    assets/Terrain/Tiles/grassDark_DM.dds
    assets/Terrain/Tiles/gravel_DM.dds
    assets/Terrain/Tiles/gravel_NM.dds
    assets/Terrain/Tiles/Snow_NM.jpg
)

set(ALL_SHADERS ${TERRAIN_SHADERS} ${SHADERS})
set(ASSETS ${TERRAIN_TEXTURES} ${TERRAIN_TILES})




add_sample_app("Atmosphere" "DiligentSamples/Samples" "${SOURCE}" "${INCLUDE}" "${ALL_SHADERS}" "${ASSETS}")

set_source_files_properties(${TERRAIN_SHADERS} PROPERTIES
    VS_DEPLOYMENT_LOCATION "shaders/terrain"
    MACOSX_PACKAGE_LOCATION "Resources/shaders/terrain"
)
set_source_files_properties(${SHADERS} PROPERTIES
    VS_DEPLOYMENT_LOCATION "shaders"
    MACOSX_PACKAGE_LOCATION "Resources/shaders"
)
set_source_files_properties(${TERRAIN_TEXTURES} PROPERTIES
    VS_DEPLOYMENT_LOCATION "Terrain"
    MACOSX_PACKAGE_LOCATION "Resources/Terrain"
)
set_source_files_properties(${TERRAIN_TILES} PROPERTIES
    VS_DEPLOYMENT_LOCATION "Terrain/Tiles"
    MACOSX_PACKAGE_LOCATION "Resources/Terrain/Tiles"
)

target_include_directories(Atmosphere
PRIVATE
    src
    src/Terrain

	${DILIGENT_ENGINE_PATH}/DiligentCore/Common/interface
	${DILIGENT_ENGINE_PATH}/DiligentCore/Graphics/GraphicsEngine/interface
	${DILIGENT_ENGINE_PATH}/DiligentCore/Graphics/GraphicsEngineVulkan/interface
	${DILIGENT_ENGINE_PATH}/DiligentCore/Graphics/GraphicsTools/interface
	${DILIGENT_ENGINE_PATH}/DiligentCore/Graphics/GraphicsAccessories/interface
	${DILIGENT_ENGINE_PATH}/DiligentCore/Primitives/interface
	${DILIGENT_ENGINE_PATH}/DiligentCore/Platforms/interface
	
	${DILIGENT_ENGINE_PATH}/DiligentFX
	${DILIGENT_ENGINE_PATH}/DiligentFX/PostProcess/EpipolarLightScattering/interface
	${DILIGENT_ENGINE_PATH}/DiligentFX/Components/interface
	${DILIGENT_ENGINE_PATH}/DiligentFX/Shaders/Common/public/
    ${DILIGENT_ENGINE_PATH}/DiligentFX/Shaders/PostProcess/EpipolarLightScattering/public
    ${DILIGENT_ENGINE_PATH}/DiligentFX/Shaders/PostProcess/ToneMapping/public
	
	${DILIGENT_ENGINE_PATH}/DiligentTools/TextureLoader/interface
	${DILIGENT_ENGINE_PATH}/DiligentTools/ThirdParty/imgui
	${DILIGENT_ENGINE_PATH}/DiligentTools/ThirdParty/imGuIZMO.quat
	${DILIGENT_ENGINE_PATH}/DiligentTools/Imgui/interface

	
	${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/include
	
)





link_directories(
	${DILIGENT_NATIVEAPPBASE_PATH}
)

file(GLOB_RECURSE  DILIGENT_TOOLS_LIB_PATH ${DILIGENT_ENGINE_PATH}/build/Win64/DiligentTools/*lib)
foreach(temp ${DILIGENT_TOOLS_LIB_PATH})
	message(${temp})
endforeach()

file(GLOB_RECURSE  Diligent_FX_LIB_PATH ${DILIGENT_ENGINE_PATH}/build/Win64/DiligentFX/*lib)
foreach(temp ${Diligent_FX_LIB_PATH})
	message(${temp})
endforeach()


target_link_libraries(Atmosphere
PRIVATE
    ${DILIGENT_TOOLS_LIB_PATH}
    ${Diligent_FX_LIB_PATH}
)

add_custom_command(TARGET Atmosphere PRE_BUILD 
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../DiligentEngine/DiligentFX/Shaders/PostProcess/EpipolarLightScattering/public/EpipolarLightScatteringStructures.fxh" "${CMAKE_CURRENT_SOURCE_DIR}/assets/shaders"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../DiligentEngine/DiligentFX/Shaders/PostProcess/EpipolarLightScattering/public/EpipolarLightScatteringFunctions.fxh" "${CMAKE_CURRENT_SOURCE_DIR}/assets/shaders"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../DiligentEngine/DiligentFX/Shaders/PostProcess/ToneMapping/public/ToneMappingStructures.fxh" "${CMAKE_CURRENT_SOURCE_DIR}/assets/shaders"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../DiligentEngine/DiligentFX/Shaders/Common/public/BasicStructures.fxh" "${CMAKE_CURRENT_SOURCE_DIR}/assets/shaders"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../DiligentEngine/DiligentFX/Shaders/Common/public/Shadows.fxh" "${CMAKE_CURRENT_SOURCE_DIR}/assets/shaders"
)

source_group("assets\\shaders\\terrain" FILES ${TERRAIN_SHADERS})
source_group("assets\\shaders" FILES ${SHADERS})
source_group("assets\\Terrain" FILES ${TERRAIN_TEXTURES})
source_group("assets\\Terrain\\Tiles" FILES ${TERRAIN_TILES})

add_definitions(-D PLATFORM_WIN32)
