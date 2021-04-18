set(PLATFORM_WIN32 False CACHE INTERNAL "")
set(PLATFORM_UNIVERSAL_WINDOWS FALSE CACHE INTERNAL "")
set(PLATFORM_ANDROID FALSE CACHE INTERNAL "")
set(PLATFORM_LINUX FALSE CACHE INTERNAL "")
set(PLATFORM_MACOS FALSE CACHE INTERNAL "")
set(PLATFORM_IOS FALSE CACHE INTERNAL "")
set(D3D11_SUPPORTED FALSE CACHE INTERNAL "D3D11 is not supported")
set(D3D12_SUPPORTED FALSE CACHE INTERNAL "D3D12 is not supported")
set(GL_SUPPORTED FALSE CACHE INTERNAL "GL is not supported")
set(GLES_SUPPORTED FALSE CACHE INTERNAL "GLES is not supported")
set(VULKAN_SUPPORTED FALSE CACHE INTERNAL "Vulkan is not supported")
set(METAL_SUPPORTED FALSE CACHE INTERNAL "Metal is not supported")


message("name" ${CMAKE_SYSTEM_NAME})
if(WIN32)
    if(${CMAKE_SYSTEM_NAME} STREQUAL "WindowsStore")
        set(PLATFORM_UNIVERSAL_WINDOWS TRUE CACHE INTERNAL "Target platform: Windows Store")
        message("Target platform: Universal Windows. SDK Version: " ${CMAKE_SYSTEM_VERSION})
    else()
        set(PLATFORM_WIN32 TRUE CACHE INTERNAL "Target platform: Win32") #WIN32 is a variable, so we cannot use string "WIN32"
        message("Target platform: Win32. SDK Version: " ${CMAKE_SYSTEM_VERSION})
    endif()
else()
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Android")
        set(PLATFORM_ANDROID TRUE CACHE INTERNAL "Target platform: Android")
        message("Target platform: Android")
    elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
        set(PLATFORM_LINUX TRUE CACHE INTERNAL "Target platform: Linux")
        message("Target Platform: Linux")
    elseif(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        if(IOS)
            set(PLATFORM_IOS TRUE CACHE INTERNAL "Target platform: iOS")
            message("Target Platform: iOS")
        else()
            set(PLATFORM_MACOS TRUE CACHE INTERNAL "Target platform: MacOS")
            message("Target Platform: MacOS")
        endif()
    elseif(${CMAKE_SYSTEM_NAME} STREQUAL "iOS")
        set(PLATFORM_IOS TRUE CACHE INTERNAL "Target platform: iOS")
        message("Target Platform: iOS")
    else()
        message(FATAL_ERROR "Unsupported platform")
    endif()
endif(WIN32)

if(PLATFORM_WIN32)
    if(MSVC)
        set(D3D11_SUPPORTED TRUE CACHE INTERNAL "D3D11 is supported on Win32 platform")
        if(CMAKE_SYSTEM_VERSION VERSION_GREATER_EQUAL "10.0")
            set(D3D12_SUPPORTED TRUE CACHE INTERNAL "D3D12 is supported on Win32 platform")
        endif()
    else()
        message("Building with MinGW")
        set(MINGW_BUILD TRUE CACHE INTERNAL "Building with MinGW")
        set(D3D11_SUPPORTED FALSE CACHE INTERNAL "D3D11 requires compiling with MSVC")
        set(D3D12_SUPPORTED FALSE CACHE INTERNAL "D3D12 requires compiling with MSVC")
    endif()

    set(GL_SUPPORTED TRUE CACHE INTERNAL "OpenGL is supported on Win32 platform")
	set(VULKAN_SUPPORTED TRUE CACHE INTERNAL "Vulkan is supported on Win32 platform")
elseif(PLATFORM_UNIVERSAL_WINDOWS)
    set(D3D11_SUPPORTED TRUE CACHE INTERNAL "D3D11 is supported on Universal Windows platform")
	if(CMAKE_SYSTEM_VERSION VERSION_GREATER_EQUAL "10.0")
		set(D3D12_SUPPORTED TRUE CACHE INTERNAL "D3D12 is supported on Universal Windows platform")
	endif()
elseif(PLATFORM_ANDROID)
    set(GLES_SUPPORTED TRUE CACHE INTERNAL "OpenGLES is supported on Android platform")
    set(VULKAN_SUPPORTED TRUE CACHE INTERNAL "Vulkan is supported on Android platform")
elseif(PLATFORM_LINUX)
    set(GL_SUPPORTED TRUE CACHE INTERNAL "OpenGL is supported on Linux platform")
    if(${ARCH} EQUAL 64)
        set(VULKAN_SUPPORTED TRUE CACHE INTERNAL "Vulkan is supported on Linux64 platform")
    endif()
elseif(PLATFORM_MACOS)
    set(GL_SUPPORTED TRUE CACHE INTERNAL "OpenGL is supported on MacOS platform")
    if(${DILIGENT_CORE_PRO_EXISTS})
        set(METAL_SUPPORTED TRUE CACHE INTERNAL "Metal is supported on MacOS platform")
    else()
        message("DiligentCorePro module is not found. Metal backend will be disabled")
    endif()
    set(VULKAN_SUPPORTED TRUE CACHE INTERNAL "Vulkan is enabled through MoltenVK on MacOS platform")
elseif(PLATFORM_IOS)
    set(GLES_SUPPORTED TRUE CACHE INTERNAL "OpenGLES is supported on iOS platform")
    if(${DILIGENT_CORE_PRO_EXISTS})
        set(METAL_SUPPORTED TRUE CACHE INTERNAL "Metal is supported on iOS platform")
    else()
        message("DiligentCorePro module is not found. Metal backend will be disabled")
    endif()
    if(VULKAN_SDK OR MoltenVK_FRAMEWORK)
        if(NOT MoltenVK_FRAMEWORK)
            set(MoltenVK_FRAMEWORK "${VULKAN_SDK}/MoltenVK/MoltenVK.xcframework" CACHE PATH "MoltenVK framework")
        endif()

        if(EXISTS ${MoltenVK_FRAMEWORK})
            set(VULKAN_SUPPORTED TRUE CACHE INTERNAL "Vulkan is enabled through MoltenVK on iOS platform")
        else()
            message(WARNING "${MoltenVK_FRAMEWORK} does not exist. Vulkan backend will be disabled.")
            unset(MoltenVK_FRAMEWORK CACHE)
        endif()
    else()
        message("VULKAN_SDK or MoltenVK_FRAMEWORK is not defined. Vulkan backend will be disabled.")
    endif()
else()
    message(FATAL_ERROR "No PLATFORM_XXX variable defined. Make sure that 'DiligentCore' folder is processed first")
endif()

if(PLATFORM_WIN32 OR PLATFORM_LINUX OR PLATFORM_MACOS)
    option(DILIGENT_BUILD_TESTS "Build Diligent Engine tests" OFF)
else()
    if(DILIGENT_BUILD_TESTS)
        message("Unit tests are not supported on this platform and will be disabled")
    endif()
    set(DILIGENT_BUILD_TESTS FALSE CACHE INTERNAL "Tests are not available on this platform" FORCE)
endif()


if(PLATFORM_WIN32)
    set(SOURCE 
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/src/Win32/SampleAppWin32.cpp
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/src/Win32/InputControllerWin32.cpp
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/src/SampleApp.cpp
    )
    set(INCLUDE 
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/include/SampleApp.hpp
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/include/Win32/InputControllerWin32.hpp
    )
    set(WIN32_RESOURCES
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/src/Win32/resources/Win32AppResource.h
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/src/Win32/resources/directx11-logo.bmp
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/src/Win32/resources/directx12-logo.bmp
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/src/Win32/resources/vulkan-logo.bmp
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/src/Win32/resources/opengl-logo.bmp
        ${DILIGENT_ENGINE_PATH}/DiligentSamples/SampleBase/src/Win32/resources/Win32AppResource.rc
    )

    function(append_sample_base_win32_source TARGET_NAME)
        # We need to add .res file to executable target to make resources available
        set(RES_FILE ${CMAKE_CURRENT_SOURCE_DIR}/../DiligentEngine/DiligentSamples/SampleBase/src/Win32/resources/Win32AppResource.rc)
        target_sources(${TARGET_NAME} PRIVATE ${RES_FILE})
        source_group("resources" FILES ${RES_FILE})
    endfunction()

elseif(PLATFORM_UNIVERSAL_WINDOWS)

    # Windows Runtime types cannot be included into static libraries
    # https://social.msdn.microsoft.com/Forums/en-US/269db513-64ef-4817-a025-43954f614eb3/lnk4264-why-are-static-libraries-not-recommended-when-authoring-windows-runtime-types?forum=winappswithnativecode
    # So as a workaround, we will include all source files into the target app project
    function(append_sample_base_uwp_source TARGET_NAME)
        get_target_property(SAMPLE_BASE_SOURCE_DIR Diligent-SampleBase SOURCE_DIR)

        set(SAMPLE_BASE_UWP_SOURCE
            ${SAMPLE_BASE_SOURCE_DIR}/src/UWP/ImguiUWPEventHelper.cpp
            ${SAMPLE_BASE_SOURCE_DIR}/src/UWP/SampleAppUWP.cpp
            ${SAMPLE_BASE_SOURCE_DIR}/src/UWP/InputControllerEventHandlerUWP.cpp
            ${SAMPLE_BASE_SOURCE_DIR}/src/SampleApp.cpp
        )

        set(SAMPLE_BASE_UWP_INCLUDE
            ${SAMPLE_BASE_SOURCE_DIR}/src/UWP/ImguiUWPEventHelper.h
            ${SAMPLE_BASE_SOURCE_DIR}/src/UWP/InputControllerEventHandlerUWP.h
            ${SAMPLE_BASE_SOURCE_DIR}/include/SampleApp.hpp
            ${SAMPLE_BASE_SOURCE_DIR}/include/UWP/InputControllerUWP.hpp
        )

        set(SAMPLE_BASE_UWP_INCLUDE_DIR
            ${SAMPLE_BASE_SOURCE_DIR}/src/UWP
        )

        target_sources(${TARGET_NAME} PRIVATE ${SAMPLE_BASE_UWP_SOURCE} ${SAMPLE_BASE_UWP_INCLUDE})
        source_group("src\\SampleBase" FILES ${SAMPLE_BASE_UWP_SOURCE})
        source_group("include\\SampleBase" FILES ${SAMPLE_BASE_UWP_INCLUDE})
        target_include_directories(${TARGET_NAME} PRIVATE ${SAMPLE_BASE_UWP_INCLUDE_DIR})
    endfunction()

elseif(PLATFORM_ANDROID)
    set(SOURCE
        src/Android/InputControllerAndroid.cpp
        src/Android/SampleAppAndroid.cpp
        src/SampleApp.cpp
    )
    set(INCLUDE 
        include/Android/InputControllerAndroid.hpp
        include/SampleApp.hpp
    )
elseif(PLATFORM_LINUX)
    set(SOURCE 
        src/Linux/InputControllerLinux.cpp
        src/Linux/SampleAppLinux.cpp
        src/SampleApp.cpp
    )
    set(INCLUDE 
        include/Linux/InputControllerLinux.hpp
        include/SampleApp.hpp
    )
elseif(PLATFORM_MACOS)

    set(SOURCE
        src/MacOS/InputControllerMacOS.cpp
        src/MacOS/SampleAppMacOS.mm
        src/SampleApp.cpp
    )
    set(INCLUDE
        Include/MacOS/InputControllerMacOS.hpp
        include/SampleApp.hpp
    )

elseif(PLATFORM_IOS)
    set(SOURCE
        src/IOS/InputControllerIOS.cpp
        src/IOS/SampleAppIOS.cpp
        src/SampleApp.cpp
    )
    set(INCLUDE
        include/IOS/InputControllerIOS.hpp
        include/SampleApp.hpp
    )

endif()


if(EXISTS ${DXC_SPIRV_PATH})
    set(DILIGENT_DXCOMPILER_FOR_SPIRV_PATH "${DXC_SPIRV_PATH}" CACHE INTERNAL "" FORCE)
    message(STATUS "Found DXCompiler for Vulkan")
else()
    set(DILIGENT_DXCOMPILER_FOR_SPIRV_PATH "" CACHE INTERNAL "" FORCE)
endif()

if(PLATFORM_WIN32 OR PLATFORM_UNIVERSAL_WINDOWS)

    function(copy_required_dlls TARGET_NAME)
		message("copy required dlls")
        if(D3D11_SUPPORTED)
            list(APPEND ENGINE_DLLS Diligent-GraphicsEngineD3D11-shared)
        endif()
        if(D3D12_SUPPORTED)
            list(APPEND ENGINE_DLLS Diligent-GraphicsEngineD3D12-shared)
        endif()
        if(GL_SUPPORTED)
            list(APPEND ENGINE_DLLS Diligent-GraphicsEngineOpenGL-shared)
        endif()
        if(VULKAN_SUPPORTED)
            list(APPEND ENGINE_DLLS Diligent-GraphicsEngineVk-shared)
        endif()
        if(METAL_SUPPORTED)
            list(APPEND ENGINE_DLLS Diligent-GraphicsEngineMetal-shared)
        endif()

		string(REGEX REPLACE ".*/\(.*\)" "\\1" CURDIR ${CMAKE_CURRENT_SOURCE_DIR})
		message("DIR " ${CMAKE_CURRENT_SOURCE_DIR})
		
		message(${TARGET_NAME})
        foreach(DLL ${ENGINE_DLLS})
			message("DLL " ${DLL})
            add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    "\"$<TARGET_FILE:${DLL}>\""
                    "\"$<TARGET_FILE_DIR:${TARGET_NAME}>\"")
        endforeach(DLL)

        # Copy D3Dcompiler_47.dll, dxcompiler.dll, and dxil.dll
        if(MSVC)
            if ((D3D11_SUPPORTED OR D3D12_SUPPORTED) AND VS_D3D_COMPILER_PATH)
                # Note that VS_D3D_COMPILER_PATH can only be used in a Visual Studio command
                # and is not a valid path during CMake configuration
                list(APPEND SHADER_COMPILER_DLLS ${VS_D3D_COMPILER_PATH})
            endif()

            if(D3D12_SUPPORTED AND VS_DXC_COMPILER_PATH AND VS_DXIL_SIGNER_PATH)
                # For the compiler to sign the bytecode, you have to have a copy of dxil.dll in 
                # the same folder as the dxcompiler.dll at runtime.

                # Note that VS_DXC_COMPILER_PATH and VS_DXIL_SIGNER_PATH can only be used in a Visual Studio command
                # and are not valid paths during CMake configuration
                list(APPEND SHADER_COMPILER_DLLS ${VS_DXC_COMPILER_PATH})
                list(APPEND SHADER_COMPILER_DLLS ${VS_DXIL_SIGNER_PATH})
            endif()

            foreach(DLL ${SHADER_COMPILER_DLLS})
                add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_if_different
                        ${DLL}
                        "\"$<TARGET_FILE_DIR:${TARGET_NAME}>\"")
            endforeach(DLL)

            if(VULKAN_SUPPORTED)
                if(NOT DEFINED DILIGENT_DXCOMPILER_FOR_SPIRV_PATH)
                    message(FATAL_ERROR "DILIGENT_DXCOMPILER_FOR_SPIRV_PATH is undefined, check order of cmake includes")
                endif()
                if(EXISTS ${DILIGENT_DXCOMPILER_FOR_SPIRV_PATH})
                    add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
                        COMMAND ${CMAKE_COMMAND} -E copy_if_different
                            ${DILIGENT_DXCOMPILER_FOR_SPIRV_PATH}
                            "\"$<TARGET_FILE_DIR:${TARGET_NAME}>/spv_dxcompiler.dll\"")
                endif()
            endif()
        endif()
    endfunction()

    function(package_required_dlls TARGET_NAME)
        if(D3D12_SUPPORTED AND VS_DXC_COMPILER_PATH AND VS_DXIL_SIGNER_PATH)
            # Copy the dlls to the project's CMake binary dir

            # Note that VS_DXC_COMPILER_PATH and VS_DXIL_SIGNER_PATH can only be used in a Visual Studio command
            # and are not valid paths during CMake configuration
            add_custom_command(TARGET ${TARGET_NAME} PRE_BUILD
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${VS_DXC_COMPILER_PATH}
                    "\"${CMAKE_CURRENT_BINARY_DIR}/dxcompiler.dll\""
                COMMAND ${CMAKE_COMMAND} -E copy_if_different
                    ${VS_DXIL_SIGNER_PATH}
                    "\"${CMAKE_CURRENT_BINARY_DIR}/dxil.dll\"")
            set(DLLS "${CMAKE_CURRENT_BINARY_DIR}/dxcompiler.dll" "${CMAKE_CURRENT_BINARY_DIR}/dxil.dll")

            # Add the dlls to the target project as source files
            target_sources(${TARGET_NAME} PRIVATE ${DLLS})

            # Label them as content
            set_source_files_properties(${DLLS} PROPERTIES 
                GENERATED TRUE
                VS_DEPLOYMENT_CONTENT 1
                VS_DEPLOYMENT_LOCATION ".")
        endif()
    endfunction()

    # Set dll output name by adding _{32|64}{r|d} suffix
    function(set_dll_output_name TARGET_NAME OUTPUT_NAME_WITHOUT_SUFFIX)
        foreach(DBG_CONFIG ${DEBUG_CONFIGURATIONS})
            set_target_properties(${TARGET_NAME} PROPERTIES
                OUTPUT_NAME_${DBG_CONFIG} ${OUTPUT_NAME_WITHOUT_SUFFIX}${DLL_DBG_SUFFIX}
            )
        endforeach()

        foreach(REL_CONFIG ${RELEASE_CONFIGURATIONS})
            set_target_properties(${TARGET_NAME} PROPERTIES
                OUTPUT_NAME_${REL_CONFIG} ${OUTPUT_NAME_WITHOUT_SUFFIX}${DLL_REL_SUFFIX}
            )
        endforeach()
    endfunction()

endif(PLATFORM_WIN32 OR PLATFORM_UNIVERSAL_WINDOWS)


function(set_common_target_properties TARGET)

    if(COMMAND custom_pre_configure_target)
        custom_pre_configure_target(${TARGET})
        if(TARGET_CONFIGURATION_COMPLETE)
            return()
        endif()
    endif()

    get_target_property(TARGET_TYPE ${TARGET} TYPE)

    if(MSVC)
        # For msvc, enable link-time code generation for release builds (I was not able to 
        # find any way to set these settings through interface library BuildSettings)
        if(TARGET_TYPE STREQUAL STATIC_LIBRARY)

            foreach(REL_CONFIG ${RELEASE_CONFIGURATIONS})
                set_target_properties(${TARGET} PROPERTIES
                    STATIC_LIBRARY_FLAGS_${REL_CONFIG} /LTCG
                )
            endforeach()

        else()

            foreach(REL_CONFIG ${RELEASE_CONFIGURATIONS})
                set_target_properties(${TARGET} PROPERTIES
                    LINK_FLAGS_${REL_CONFIG} "/LTCG /OPT:REF /INCREMENTAL:NO"
                )
            endforeach()

            if(PLATFORM_UNIVERSAL_WINDOWS)
                # On UWP, disable incremental link to avoid linker warnings
                foreach(DBG_CONFIG ${DEBUG_CONFIGURATIONS})
                    set_target_properties(${TARGET} PROPERTIES
                        LINK_FLAGS_${DBG_CONFIG} "/INCREMENTAL:NO"
                    )
                endforeach()
            endif()
        endif()
    else()
        set_target_properties(${TARGET} PROPERTIES
            CXX_VISIBILITY_PRESET hidden # -fvisibility=hidden
            C_VISIBILITY_PRESET hidden # -fvisibility=hidden
            VISIBILITY_INLINES_HIDDEN TRUE

            # Without -fPIC option GCC fails to link static libraries into dynamic library:
            #  -fPIC  
            #      If supported for the target machine, emit position-independent code, suitable for 
            #      dynamic linking and avoiding any limit on the size of the global offset table.
            POSITION_INDEPENDENT_CODE ON

            # It is crucial to set CXX_STANDARD flag to only affect c++ files and avoid failures compiling c-files:
            # error: invalid argument '-std=c++11' not allowed with 'C/ObjC'
            CXX_STANDARD 11
            CXX_STANDARD_REQUIRED ON

            C_STANDARD 11
        )

        if(NOT MINGW_BUILD)
            # Do not disable extensions when building with MinGW!
            set_target_properties(${TARGET} PROPERTIES
                CXX_EXTENSIONS OFF
            )
        endif()
    endif()

    if(COMMAND custom_post_configure_target)
        custom_post_configure_target(${TARGET})
    endif()

endfunction()

function(find_targets_in_directory _RESULT _DIR)
    get_property(_subdirs DIRECTORY "${_DIR}" PROPERTY SUBDIRECTORIES)
    foreach(_subdir IN LISTS _subdirs)
        find_targets_in_directory(${_RESULT} "${_subdir}")
    endforeach()
    get_property(_SUB_TARGETS DIRECTORY "${_DIR}" PROPERTY BUILDSYSTEM_TARGETS)
    set(${_RESULT} ${${_RESULT}} ${_SUB_TARGETS} PARENT_SCOPE)
endfunction()

function(set_directory_root_folder _DIRECTORY _ROOT_FOLDER)
    find_targets_in_directory(_TARGETS ${_DIRECTORY})
    foreach(_TARGET IN LISTS _TARGETS)
        get_target_property(_FOLDER ${_TARGET} FOLDER)
        if(_FOLDER)
            set_target_properties(${_TARGET} PROPERTIES FOLDER "${_ROOT_FOLDER}/${_FOLDER}")
        else()
            set_target_properties(${_TARGET} PROPERTIES FOLDER "${_ROOT_FOLDER}")
        endif()
    endforeach()
endfunction()


# Returns default backend library type (static/dynamic) for the current platform
function(get_backend_libraries_type _LIB_TYPE)
    if(PLATFORM_WIN32 OR PLATFORM_LINUX OR PLATFORM_ANDROID OR PLATFORM_UNIVERSAL_WINDOWS OR PLATFORM_MACOS)
        set(LIB_TYPE "shared")
    elseif(PLATFORM_IOS)
        # Statically link with the engine on iOS.
        # It is also possible to link dynamically by
        # putting the library into the framework.
        set(LIB_TYPE "static")
    else()
        message(FATAL_ERROR "Undefined platform")
    endif()
    set(${_LIB_TYPE} ${LIB_TYPE} PARENT_SCOPE)
endfunction()


# Adds the list of supported backend targets to variable ${_TARGETS} in parent scope.
# Second argument to the function may override the target type (static/dynamic). If It
# is not given, default target type for the platform is used.
function(get_supported_backends _TARGETS)
    if(${ARGC} GREATER 1)
        set(LIB_TYPE ${ARGV1})
    else()
        get_backend_libraries_type(LIB_TYPE)
    endif()

    if(D3D11_SUPPORTED)
	    list(APPEND BACKENDS Diligent-GraphicsEngineD3D11-${LIB_TYPE})
    endif()
    if(D3D12_SUPPORTED)
	    list(APPEND BACKENDS Diligent-GraphicsEngineD3D12-${LIB_TYPE})
    endif()
    if(GL_SUPPORTED OR GLES_SUPPORTED)
	    list(APPEND BACKENDS Diligent-GraphicsEngineOpenGL-${LIB_TYPE})
    endif()
    if(VULKAN_SUPPORTED)
	    list(APPEND BACKENDS Diligent-GraphicsEngineVk-${LIB_TYPE})
    endif()
    if(METAL_SUPPORTED)
	    list(APPEND BACKENDS Diligent-GraphicsEngineMetal-${LIB_TYPE})
    endif()
    # ${_TARGETS} == ENGINE_LIBRARIES
    # ${${_TARGETS}} == ${ENGINE_LIBRARIES}
    set(${_TARGETS} ${${_TARGETS}} ${BACKENDS} PARENT_SCOPE)
endfunction()


# Returns path to the target relative to CMake root
function(get_target_relative_dir _TARGET _DIR)
    get_target_property(TARGET_SOURCE_DIR ${_TARGET} SOURCE_DIR)
    file(RELATIVE_PATH TARGET_RELATIVE_PATH "${CMAKE_SOURCE_DIR}" "${TARGET_SOURCE_DIR}")
    set(${_DIR} ${TARGET_RELATIVE_PATH} PARENT_SCOPE)
endfunction()

# Performs installation steps for the core library
function(install_core_lib _TARGET)
    get_target_relative_dir(${_TARGET} TARGET_RELATIVE_PATH)

    get_target_property(TARGET_TYPE ${_TARGET} TYPE)
    if(TARGET_TYPE STREQUAL STATIC_LIBRARY)
        list(APPEND DILIGENT_CORE_INSTALL_LIBS_LIST ${_TARGET})
        set(DILIGENT_CORE_INSTALL_LIBS_LIST ${DILIGENT_CORE_INSTALL_LIBS_LIST} CACHE INTERNAL "Core libraries installation list")
    elseif(TARGET_TYPE STREQUAL SHARED_LIBRARY)
        install(TARGETS				 ${_TARGET}
                ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}/${DILIGENT_CORE_DIR}/$<CONFIG>"
                LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}/${DILIGENT_CORE_DIR}/$<CONFIG>"
                RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}/${DILIGENT_CORE_DIR}/$<CONFIG>"
        )
        if (DILIGENT_INSTALL_PDB)
            install(FILES $<TARGET_PDB_FILE:${_TARGET}> DESTINATION "${CMAKE_INSTALL_BINDIR}/${DILIGENT_CORE_DIR}/$<CONFIG>" OPTIONAL)
        endif()
    endif()

    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/interface")
        install(DIRECTORY    interface
                DESTINATION  "${CMAKE_INSTALL_INCLUDEDIR}/${TARGET_RELATIVE_PATH}/"
        )
    endif()
endfunction()


function(install_combined_static_lib COMBINED_LIB_NAME LIBS_LIST CUSTOM_TARGET_NAME CUSTOM_TARGET_FOLDER INSTALL_DESTINATION)

    foreach(LIB ${LIBS_LIST})
        list(APPEND COMBINED_LIB_TARGET_FILES $<TARGET_FILE:${LIB}>)
    endforeach(LIB)

    if(MSVC)
        add_custom_command(
            OUTPUT ${COMBINED_LIB_NAME}
            COMMAND lib.exe /OUT:${COMBINED_LIB_NAME} ${COMBINED_LIB_TARGET_FILES}
            DEPENDS ${LIBS_LIST}
            COMMENT "Combining libraries..."
        )
        add_custom_target(${CUSTOM_TARGET_NAME} ALL DEPENDS ${COMBINED_LIB_NAME})
    else()

        if(PLATFORM_WIN32)
            # do NOT use stock ar on MinGW
            find_program(AR NAMES x86_64-w64-mingw32-gcc-ar)
        else()
            set(AR ${CMAKE_AR})
        endif()

        if(AR)
            add_custom_command(
                OUTPUT ${COMBINED_LIB_NAME}
                # Delete all object files from current directory
                COMMAND ${CMAKE_COMMAND} -E remove "*${CMAKE_C_OUTPUT_EXTENSION}"
                DEPENDS ${LIBS_LIST}
                COMMENT "Combining libraries..."
            )

            # Unpack all object files from all targets to current directory
            foreach(LIB_TARGET ${COMBINED_LIB_TARGET_FILES})
                add_custom_command(
                    OUTPUT ${COMBINED_LIB_NAME}
                    COMMAND ${AR} -x ${LIB_TARGET}
                    APPEND
                )
            endforeach()

            # Pack object files to a combined library and delete them
            add_custom_command(
                OUTPUT ${COMBINED_LIB_NAME}
                COMMAND ${AR} -crs ${COMBINED_LIB_NAME} "*${CMAKE_C_OUTPUT_EXTENSION}"
                COMMAND ${CMAKE_COMMAND} -E remove "*${CMAKE_C_OUTPUT_EXTENSION}"
                APPEND
            )

            add_custom_target(${CUSTOM_TARGET_NAME} ALL DEPENDS ${COMBINED_LIB_NAME})
        else()
            message("ar command is not found")
        endif()
    endif()

    if(TARGET ${CUSTOM_TARGET_NAME})
        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${COMBINED_LIB_NAME}"
                DESTINATION ${INSTALL_DESTINATION}
        )
        set_target_properties(${CUSTOM_TARGET_NAME} PROPERTIES
            FOLDER ${CUSTOM_TARGET_FOLDER}
        )
    else()
        message("Unable to find librarian tool. Combined ${COMBINED_LIB_NAME} static library will not be produced.")
    endif()

endfunction()




function(add_format_validation_target MODULE_NAME MODULE_ROOT_PATH IDE_FOLDER)

    if(${DILIGENT_NO_FORMAT_VALIDATION})
        return()
    endif()

    # Start by copying .clang-format file to the module's root folder
    add_custom_target(${MODULE_NAME}-ValidateFormatting ALL
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${DILIGENT_CORE_SOURCE_DIR}/.clang-format" "${MODULE_ROOT_PATH}/.clang-format"
    )

    if(CMAKE_HOST_SYSTEM_NAME STREQUAL "Windows")
        set(RUN_VALIDATION_SCRIPT validate_format_win.bat)
    elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Linux")
        set(RUN_VALIDATION_SCRIPT ./validate_format_linux.sh)
    elseif(CMAKE_HOST_SYSTEM_NAME STREQUAL "Darwin")
        set(RUN_VALIDATION_SCRIPT ./validate_format_mac.sh)
    else()
        mesage(FATAL_ERROR "Unexpected host system")
    endif()

    # Run the format validation script
    add_custom_command(TARGET ${MODULE_NAME}-ValidateFormatting
        COMMAND ${RUN_VALIDATION_SCRIPT}
        WORKING_DIRECTORY "${MODULE_ROOT_PATH}/BuildTools/FormatValidation"
        COMMENT "Validating ${MODULE_NAME} module's source code formatting..."
        VERBATIM
    )

    if(TARGET ${MODULE_NAME}-ValidateFormatting)
        set_target_properties(${MODULE_NAME}-ValidateFormatting PROPERTIES FOLDER ${IDE_FOLDER})
    endif()

endfunction()



if(PLATFORM_WIN32)

    set(SOURCE 
        ${DILIGENT_ENGINE_PATH}/DiligentTools/NativeApp/src/Win32/WinMain.cpp
    )
    set(INCLUDE
        ${DILIGENT_ENGINE_PATH}/DiligentTools/NativeApp/include/Win32/Win32AppBase.hpp
    )

    function(add_win32_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_executable(${TARGET_NAME} WIN32 ${SOURCE} ${INCLUDE} ${ASSETS})

        if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
            # libmingw32 must be included BEFORE Diligent-NativeAppBase that contains the definition of WinMain.
            # otherwise WinMain will be stripped out of Diligent-NativeAppBase and will be unresolved.
            target_link_libraries(${TARGET_NAME}
            PRIVATE
                mingw32
            )
        endif()
    endfunction()

    function(add_target_platform_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_win32_app("${TARGET_NAME}" "${SOURCE}" "${INCLUDE}" "${ASSETS}")
    endfunction()

elseif(PLATFORM_UNIVERSAL_WINDOWS)

    set(SOURCE 
        src/UWP/dummy.cpp
    )
    set(INCLUDE
        include/UWP/UWPAppBase.hpp
    )

    # Windows Runtime types cannot be included into static libraries
    # https://social.msdn.microsoft.com/Forums/en-US/269db513-64ef-4817-a025-43954f614eb3/lnk4264-why-are-static-libraries-not-recommended-when-authoring-windows-runtime-types?forum=winappswithnativecode
    # So as a workaround, we will include all source files into the target app project
    function(add_uwp_app TARGET_NAME SOURCE INCLUDE ASSETS)
        get_target_property(NATIVE_APP_SOURCE_DIR Diligent-NativeAppBase SOURCE_DIR)

        set(UWP_SOURCE
            ${NATIVE_APP_SOURCE_DIR}/src/UWP/Common/DeviceResources.cpp
            ${NATIVE_APP_SOURCE_DIR}/src/UWP/App.cpp
            ${NATIVE_APP_SOURCE_DIR}/src/UWP/UWPAppBase.cpp
        )

        set(UWP_INCLUDE
            ${NATIVE_APP_SOURCE_DIR}/src/UWP/Common/DeviceResources.h
            ${NATIVE_APP_SOURCE_DIR}/src/UWP/Common/DirectXHelper.h
            ${NATIVE_APP_SOURCE_DIR}/src/UWP/App.h
            ${NATIVE_APP_SOURCE_DIR}/include/UWP/UWPAppBase.hpp
            ${NATIVE_APP_SOURCE_DIR}/src/UWP/Common/StepTimer.h
        )

        add_executable(${TARGET_NAME} WIN32 ${SOURCE} ${INCLUDE} ${ASSETS} ${UWP_SOURCE} ${UWP_INCLUDE})
        set_source_files_properties(${ASSETS} PROPERTIES VS_DEPLOYMENT_CONTENT 1)

        target_include_directories(${TARGET_NAME} 
        PUBLIC
            ${NATIVE_APP_SOURCE_DIR}/Src/UWP
        )
    
        target_link_libraries(${TARGET_NAME} PRIVATE dxgi.lib)

        source_group("UWP Common\\src" FILES ${UWP_SOURCE})
        source_group("UWP Common\\include" FILES ${UWP_INCLUDE})

    endfunction(add_uwp_app)

    function(add_target_platform_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_uwp_app("${TARGET_NAME}" "${SOURCE}" "${INCLUDE}" "${ASSETS}")
    endfunction()


elseif(PLATFORM_ANDROID)

    set(SOURCE 
        src/Android/AndroidAppBase.cpp
    )
    set(INCLUDE
        include/Android/AndroidAppBase.hpp
    )
    function(add_android_app TARGET_NAME SOURCE INCLUDE ASSETS)
        get_target_property(NATIVE_APP_SOURCE_DIR Diligent-NativeAppBase SOURCE_DIR)
        set(ANDROID_SOURCE
            ${NATIVE_APP_SOURCE_DIR}/src/Android/AndroidMain.cpp
        )
        add_library(${TARGET_NAME} SHARED ${SOURCE} ${ANDROID_SOURCE} ${INCLUDE} ${ASSETS})
        target_link_libraries(${TARGET_NAME} 
        PRIVATE 
            android
            native_app_glue
        )
        # Export ANativeActivity_onCreate(),
        # Refer to: https://github.com/android-ndk/ndk/issues/381.
        set_target_properties(${TARGET_NAME}
        PROPERTIES
          LINK_FLAGS "-u ANativeActivity_onCreate"
        )
        #target_include_directories(${TARGET_NAME} 
        #PRIVATE
        #	${ANDROID_NDK}/sources/android/cpufeatures
        #)
        source_group("Android" FILES ${ANDROID_SOURCE})
    endfunction()

    function(add_target_platform_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_android_app("${TARGET_NAME}" "${SOURCE}" "${INCLUDE}" "${ASSETS}")
    endfunction()

elseif(PLATFORM_LINUX)

    set(SOURCE 
        src/Linux/LinuxMain.cpp
    )
    set(INCLUDE
        include/Linux/LinuxAppBase.hpp
    )
    function(add_linux_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_executable(${TARGET_NAME} ${SOURCE} ${INCLUDE} ${ASSETS})
    endfunction()

    function(add_target_platform_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_linux_app("${TARGET_NAME}" "${SOURCE}" "${INCLUDE}" "${ASSETS}")
    endfunction()


elseif(PLATFORM_MACOS)

    set(SOURCE
        src/MacOS/MacOSAppBase.cpp
    )
    set(INCLUDE
        include/MacOS/MacOSAppBase.hpp
    )

    function(add_macos_app TARGET_NAME SOURCE INCLUDE ASSETS)
        get_target_property(NATIVE_APP_SOURCE_DIR Diligent-NativeAppBase SOURCE_DIR)

        set(APPLE_SOURCE
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/main.m
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/WindowController.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/AppDelegate.m
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/FullscreenWindow.m
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/GLView.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/MetalView.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/MVKView.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/ViewBase.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/ViewController.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/ModeSelectionViewController.mm
        )

        set(APPLE_INCLUDE
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/WindowController.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/AppDelegate.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/FullscreenWindow.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/GLView.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/MetalView.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/MVKView.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/ViewBase.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/ViewController.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX/ModeSelectionViewController.h
        )

        set(APPLE_RESOURCES
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/OSX/Base.lproj/Main.storyboard
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/OSX/Images.xcassets/AppIcon.appiconset/dg.icns
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/OSX/Images.xcassets/opengl-logo.png
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/OSX/Images.xcassets/vulkan-logo.png
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/OSX/Images.xcassets/metal-logo.png
        )

        set(APPLE_INFO_PLIST
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/OSX/Info.plist
        )

        set(APPLE_INCLUDE_DIRS
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/OSX
        )

        add_executable(${TARGET_NAME} MACOSX_BUNDLE ${SOURCE} ${APPLE_SOURCE} ${INCLUDE} ${APPLE_INCLUDE} ${ASSETS} ${APPLE_RESOURCES})
        set_target_properties(${TARGET_NAME} PROPERTIES
            MACOSX_BUNDLE_GUI_IDENTIFIER "com.diligentengine.samples.${TARGET_NAME}"
            MACOSX_BUNDLE_INFO_PLIST "${APPLE_INFO_PLIST}"
            RESOURCE "${APPLE_RESOURCES}"
        )
        source_group("MacOS" FILES ${APPLE_SOURCE})
        source_group("MacOS" FILES ${APPLE_INCLUDE})
        source_group("Resources" FILES ${APPLE_RESOURCES})
        target_include_directories(${TARGET_NAME} PRIVATE ${APPLE_INCLUDE_DIRS})

        find_package(OpenGL REQUIRED)

        find_library(CORE_VIDEO CoreVideo)
        if (NOT CORE_VIDEO)
                message(FATAL_ERROR "CoreVideo is not found")
        endif()

        find_library(METAL_FRAMEWORK Metal)
        if (NOT METAL_FRAMEWORK)
            message(FATAL_ERROR "Metal framework is not found")
        endif()

        find_library(CORE_ANIMATION QuartzCore)
        if (NOT CORE_ANIMATION)
            message(FATAL_ERROR "QuartzCore (CoreAnimation) is not found")
        endif()

        target_link_libraries(${TARGET_NAME} PRIVATE ${OPENGL_LIBRARY} ${CORE_VIDEO} ${METAL_FRAMEWORK} ${CORE_ANIMATION})

    endfunction()

    function(add_target_platform_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_macos_app("${TARGET_NAME}" "${SOURCE}" "${INCLUDE}" "${ASSETS}")
    endfunction()

elseif(PLATFORM_IOS)

    set(SOURCE
        src/IOS/IOSAppBase.cpp
    )
    set(INCLUDE
        include/IOS/IOSAppBase.hpp
    )

    function(add_ios_app TARGET_NAME SOURCE INCLUDE ASSETS)
        get_target_property(NATIVE_APP_SOURCE_DIR Diligent-NativeAppBase SOURCE_DIR)

        set(APPLE_SOURCE
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/main.m
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/AppDelegate.m
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/BaseView.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/EAGLView.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/AppViewBase.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/ModeSelectionViewController.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/MetalView.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/MVKView.mm
        )

        set(APPLE_INCLUDE
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/AppDelegate.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/BaseView.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/EAGLView.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/AppViewBase.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/ModeSelectionViewController.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/MetalView.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS/MVKView.h
        )

        set(APPLE_RESOURCES
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/iOS/Base.lproj/Main.storyboard
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/iOS/Base.lproj/LaunchScreen.xib
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/iOS/Images.xcassets/AppIcon.appiconset/dg-icon.png
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/iOS/Images.xcassets/opengles-logo.png
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/iOS/Images.xcassets/vulkan-logo.png
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/iOS/Images.xcassets/metal-logo.png
        )

        set(APPLE_INFO_PLIST
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/iOS/info.plist
        )

        set(APPLE_INCLUDE_DIRS
                ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/iOS
        )

        add_executable(${TARGET_NAME} MACOSX_BUNDLE ${SOURCE} ${APPLE_SOURCE} ${INCLUDE} ${APPLE_INCLUDE} ${ASSETS} ${APPLE_RESOURCES})
        set_target_properties(${TARGET_NAME} PROPERTIES
            MACOSX_BUNDLE_GUI_IDENTIFIER "com.diligentengine.samples.${TARGET_NAME}"
            MACOSX_BUNDLE_INFO_PLIST "${APPLE_INFO_PLIST}"
            RESOURCE "${APPLE_RESOURCES}"
            XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY "iPhone Developer"
            # XCODE_ATTRIBUTE_DEVELOPMENT_TEAM "Dev Team"
            BUILD_RPATH "@executable_path"
        )
        source_group("iOS" FILES ${APPLE_SOURCE})
        source_group("iOS" FILES ${APPLE_INCLUDE})
        source_group("Resources" FILES ${APPLE_RESOURCES})
        target_include_directories(${TARGET_NAME} PRIVATE ${APPLE_INCLUDE_DIRS})

        find_library(OPENGLES OpenGLES)
        if (NOT OPENGLES)
            message(FATAL_ERROR "OpenGLES is not found")
        endif()

        find_library(UIKIT UIKit)
        if (NOT UIKIT)
            message(FATAL_ERROR "UIKIT is not found")
        endif()

        find_library(CORE_ANIMATION QuartzCore)
        if (NOT CORE_ANIMATION)
            message(FATAL_ERROR "QuartzCore (CoreAnimation) is not found")
        endif()

        target_link_libraries(${TARGET_NAME} PRIVATE ${OPENGLES} ${UIKIT} ${CORE_ANIMATION})
    endfunction()

    function(add_target_platform_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_ios_app("${TARGET_NAME}" "${SOURCE}" "${INCLUDE}" "${ASSETS}")
    endfunction()

else()
    message(FATAL_ERROR "Unknown platform")
endif()


function(add_sample_app APP_NAME IDE_FOLDER SOURCE INCLUDE SHADERS ASSETS)

    set_source_files_properties(${SHADERS} PROPERTIES VS_TOOL_OVERRIDE "None")
    set(ALL_ASSETS ${ASSETS} ${SHADERS})
    add_target_platform_app(${APP_NAME} "${SOURCE}" "${INCLUDE}" "${ALL_ASSETS}")

    set_source_files_properties(${ALL_ASSETS} PROPERTIES 
        VS_DEPLOYMENT_LOCATION "."
        MACOSX_PACKAGE_LOCATION "Resources" 
    )

    if(PLATFORM_WIN32)
        set_target_properties(${APP_NAME} PROPERTIES 
            VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/assets"
        )
        copy_required_dlls(${APP_NAME})
        append_sample_base_win32_source(${APP_NAME})
    elseif(PLATFORM_UNIVERSAL_WINDOWS)
        append_sample_base_uwp_source(${APP_NAME})
        package_required_dlls(${APP_NAME})
    endif()

    target_include_directories(${APP_NAME}
    PRIVATE	
        src
    )

	find_library(DILIGENT_NATIVEAPPBASE_PATH Diligent-NativeAppBase ${DILIGENT_ENGINE_PATH}/build/Win64/DiligentTools/NativeApp/Debug)
	IF(NOT DILIGENT_NATIVEAPPBASE_PATH)
		MESSAGE(FATAL_ERROR "Diligent-NativeAppBase not found")
	ENDIF(NOT DILIGENT_NATIVEAPPBASE_PATH)
	MESSAGE(STATUS ${DILIGENT_NATIVEAPPBASE_PATH} " found")
	
	file(GLOB_RECURSE  DILIGENT_BUILDSETTINGS_PATH ${DILIGENT_ENGINE_PATH}/build/Win64/DiligentCore/*lib)
	foreach(temp ${DILIGENT_BUILDSETTINGS_PATH})
		message(${temp})
	endforeach()
	
	find_library(DILIGENT_SAMPLEBASE_PATH Diligent-SampleBase ${DILIGENT_ENGINE_PATH}/build/Win64/DiligentSamples/SampleBase/Debug)
	IF(NOT DILIGENT_SAMPLEBASE_PATH)
		MESSAGE(FATAL_ERROR "Diligent-NativeAppBase not found")
	ENDIF(NOT DILIGENT_SAMPLEBASE_PATH)
	MESSAGE(STATUS ${DILIGENT_SAMPLEBASE_PATH} " found")
	
	

    target_link_libraries(${APP_NAME}
    PRIVATE
        # On Linux we must have Diligent-NativeAppBase go first, otherwise the linker 
        # will fail to resolve Diligent::CreateApplication() function.
        ${DILIGENT_NATIVEAPPBASE_PATH} #Diligent-NativeAppBase
        ${DILIGENT_BUILDSETTINGS_PATH} #Diligent-BuildSettings
        ${DILIGENT_SAMPLEBASE_PATH} #Diligent-SampleBase
    )
    set_common_target_properties(${APP_NAME})

    if(MSVC)
        # Disable MSVC-specific warnings
        # - w4201: nonstandard extension used: nameless struct/union
        target_compile_options(${APP_NAME} PRIVATE /wd4201)
    endif()

    set_target_properties(${APP_NAME} PROPERTIES
        FOLDER ${IDE_FOLDER}
    )

    source_group("src" FILES ${SOURCE} ${INCLUDE})
    source_group("assets" FILES ${ALL_ASSETS})	

    target_sources(${APP_NAME} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/readme.md")
    set_source_files_properties(
        "${CMAKE_CURRENT_SOURCE_DIR}/readme.md" PROPERTIES HEADER_FILE_ONLY TRUE
    )

    if(PLATFORM_WIN32 OR PLATFORM_LINUX)
        # Copy assets to target folder
        add_custom_command(TARGET ${APP_NAME} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_directory
                "${CMAKE_CURRENT_SOURCE_DIR}/assets"
                "\"$<TARGET_FILE_DIR:${APP_NAME}>\"")
    endif()

    if(PLATFORM_MACOS AND VULKAN_LIB_PATH)
        # Configure rpath so that executables can find vulkan library
        set_target_properties(${APP_NAME} PROPERTIES
            BUILD_RPATH "${VULKAN_LIB_PATH}"
        )
    endif()

    if(DILIGENT_INSTALL_SAMPLES)
        # Install instructions
        file(RELATIVE_PATH TUTORIAL_REL_PATH "${CMAKE_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
        install(TARGETS     ${APP_NAME}
                DESTINATION "${CMAKE_INSTALL_BINDIR}/${TUTORIAL_REL_PATH}/$<CONFIG>")

        if(PLATFORM_LINUX OR PLATFORM_WIN32)
            install(DIRECTORY   "${CMAKE_CURRENT_SOURCE_DIR}/assets/"
                    DESTINATION "${CMAKE_INSTALL_BINDIR}/${TUTORIAL_REL_PATH}/$<CONFIG>")
        endif()

        if(PLATFORM_WIN32)
            get_supported_backends(BACKEND_LIBRARIES)
            install(TARGETS  ${BACKEND_LIBRARIES}
                    RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}/${TUTORIAL_REL_PATH}/$<CONFIG>"
                    LIBRARY DESTINATION "${CMAKE_INSTALL_BINDIR}/${TUTORIAL_REL_PATH}/$<CONFIG>"
                    ARCHIVE DESTINATION "${CMAKE_INSTALL_BINDIR}/${TUTORIAL_REL_PATH}/$<CONFIG>")
        endif()

        if(PLATFORM_LINUX)
            set_target_properties(${APP_NAME} PROPERTIES
                INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}/${DILIGENT_CORE_DIR}/${CMAKE_BUILD_TYPE}"
            )
        endif()
    endif()

endfunction()