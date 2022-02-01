cmake_minimum_required (VERSION 3.6)

project(main CXX)

Set(PLATFORM_WIN32 True)


add_compile_definitions(_DEBUG)
add_compile_definitions(DEBUG)
add_compile_definitions(PLATFORM_WIN32=1)
add_compile_definitions(D3D11_SUPPORTED=0)
add_compile_definitions(D3D12_SUPPORTED=0)
add_compile_definitions(GL_SUPPORTED=0)
add_compile_definitions(GLES_SUPPORTED=0)
add_compile_definitions(VULKAN_SUPPORTED=1)
add_compile_definitions(METAL_SUPPORTED=0)
add_compile_definitions(DILIGENT_DEVELOPMENT)
add_compile_definitions(DILIGENT_DEBUG)
add_compile_definitions(ENGINE_DLL=1)


set(DILIGENT_ENGINE_DIR D:/indie/DiligentEngine)
set(DILIGENT_ENGINE_BUILD_DIR ${DILIGENT_ENGINE_DIR}/build/Win64) 



include(BuildUtils.cmake)




if(PLATFORM_WIN32)

    set(SOURCE 
        src/Win32/WinMain.cpp
    )
    set(INCLUDE
        include/Win32/Win32AppBase.hpp
    )

    function(add_win32_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_executable(${TARGET_NAME} WIN32 ${SOURCE} ${INCLUDE} ${ASSETS})
		set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${TARGET_NAME})
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
    
        source_group("UWP Common//src" FILES ${UWP_SOURCE})
        source_group("UWP Common//include" FILES ${UWP_INCLUDE})

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
            XCODE_ATTRIBUTE_PRODUCT_BUNDLE_IDENTIFIER "com.diligentengine.samples.${TARGET_NAME}"
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
        # Silence OpenGL deprecation warnings
        target_compile_definitions(${TARGET_NAME} PUBLIC GL_SILENCE_DEPRECATION)

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
            XCODE_ATTRIBUTE_PRODUCT_BUNDLE_IDENTIFIER "com.diligentengine.samples.${TARGET_NAME}"
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
        # Silence OpenGLES deprecation warnings
        target_compile_definitions(${TARGET_NAME} PUBLIC GLES_SILENCE_DEPRECATION)

    endfunction()

    function(add_target_platform_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_ios_app("${TARGET_NAME}" "${SOURCE}" "${INCLUDE}" "${ASSETS}")
    endfunction()

elseif(PLATFORM_TVOS)

    set(SOURCE
        src/TVOS/TVOSAppBase.cpp
    )
    set(INCLUDE
        include/TVOS/TVOSAppBase.hpp
    )

    function(add_tvos_app TARGET_NAME SOURCE INCLUDE ASSETS)
        get_target_property(NATIVE_APP_SOURCE_DIR Diligent-NativeAppBase SOURCE_DIR)

        set(APPLE_SOURCE
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/main.m
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/tvOS/AppDelegate.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/tvOS/MainUIView.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/tvOS/MetalView.mm
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/tvOS/MetalViewController.mm
        )

        set(APPLE_INCLUDE
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/tvOS/AppDelegate.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/tvOS/MainUIView.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/tvOS/MetalView.h
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/tvOS/MetalViewController.h
        )

        set(APPLE_RESOURCES
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/tvOS/Base.lproj/Main.storyboard
        )

        set(APPLE_INFO_PLIST
            ${NATIVE_APP_SOURCE_DIR}/Apple/Data/tvOS/info.plist
        )

        set(APPLE_INCLUDE_DIRS
            ${NATIVE_APP_SOURCE_DIR}/Apple/Source/Classes/tvOS
        )

        add_executable(${TARGET_NAME} MACOSX_BUNDLE ${SOURCE} ${APPLE_SOURCE} ${INCLUDE} ${APPLE_INCLUDE} ${ASSETS} ${APPLE_RESOURCES})
        set_target_properties(${TARGET_NAME} PROPERTIES
            XCODE_ATTRIBUTE_PRODUCT_BUNDLE_IDENTIFIER "com.diligentengine.samples.${TARGET_NAME}"
            MACOSX_BUNDLE_INFO_PLIST "${APPLE_INFO_PLIST}"
            RESOURCE "${APPLE_RESOURCES}"
            XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY "iPhone Developer"
            # XCODE_ATTRIBUTE_DEVELOPMENT_TEAM "Dev Team"
            BUILD_RPATH "@executable_path"
        )
        source_group("tvOS" FILES ${APPLE_SOURCE})
        source_group("tvOS" FILES ${APPLE_INCLUDE})
        source_group("Resources" FILES ${APPLE_RESOURCES})
        target_include_directories(${TARGET_NAME} PRIVATE ${APPLE_INCLUDE_DIRS})

        find_library(UIKIT UIKit)
        if (NOT UIKIT)
            message(FATAL_ERROR "UIKIT is not found")
        endif()

        find_library(CORE_ANIMATION QuartzCore)
        if (NOT CORE_ANIMATION)
            message(FATAL_ERROR "QuartzCore (CoreAnimation) is not found")
        endif()

        target_link_libraries(${TARGET_NAME} PRIVATE ${UIKIT} ${CORE_ANIMATION})

    endfunction()

    function(add_target_platform_app TARGET_NAME SOURCE INCLUDE ASSETS)
        add_tvos_app("${TARGET_NAME}" "${SOURCE}" "${INCLUDE}" "${ASSETS}")
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
		#append_sample_base_win32_source(${APP_NAME})
    elseif(PLATFORM_UNIVERSAL_WINDOWS)
		#append_sample_base_uwp_source(${APP_NAME})
        package_required_dlls(${APP_NAME})
    endif()

    target_include_directories(${APP_NAME}
    PRIVATE	
        src
    )
    target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentTools/NativeApp/Debug/Diligent-NativeAppBase.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentSamples/SampleBase/Debug/Diligent-SampleBase.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/GraphicsAccessories/Debug/Diligent-GraphicsAccessories.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentTools/TextureLoader/Debug/Diligent-TextureLoader.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentFX/Debug/DiligentFX.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentTools/Imgui/Debug/Diligent-Imgui.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/GraphicsEngineD3D11/Debug/GraphicsEngineD3D11_64d.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/GraphicsEngineD3D12/Debug/GraphicsEngineD3D12_64d.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/GraphicsEngineOpenGL/Debug/GraphicsEngineOpenGL_64d.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/GraphicsEngineVulkan/Debug/GraphicsEngineVk_64d.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentTools/AssetLoader/Debug/Diligent-AssetLoader.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentTools/ThirdParty/lpng-1.6.17/Debug/LibPng.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentTools/ThirdParty/libtiff/Debug/LibTiff.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentTools/ThirdParty/zlib-1.2.8/Debug/ZLib.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentTools/ThirdParty/libjpeg-9a/Debug/LibJpeg.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/GraphicsTools/Debug/Diligent-GraphicsTools.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/GraphicsEngineD3DBase/Debug/Diligent-GraphicsEngineD3DBase.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/ShaderTools/Debug/Diligent-ShaderTools.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/HLSL2GLSLConverterLib/Debug/Diligent-HLSL2GLSLConverterLib.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Graphics/GraphicsEngine/Debug/Diligent-GraphicsEngine.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Common/Debug/Diligent-Common.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Platforms/Win32/Debug/Diligent-Win32Platform.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Platforms/Basic/Debug/Diligent-BasicPlatform.lib)
	target_link_libraries(${APP_NAME} PRIVATE Shlwapi.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/ThirdParty/SPIRV-Cross/Debug/spirv-cross-cored.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/ThirdParty/glslang/SPIRV/Debug/SPIRVd.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/ThirdParty/SPIRV-Tools/source/opt/Debug/SPIRV-Tools-opt.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/ThirdParty/SPIRV-Tools/source/Debug/SPIRV-Tools.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/ThirdParty/glslang/glslang/Debug/glslangd.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/ThirdParty/glslang/glslang/Debug/MachineIndependentd.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/ThirdParty/glslang/glslang/Debug/GenericCodeGend.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/ThirdParty/glslang/OGLCompilersDLL/Debug/OGLCompilerd.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/ThirdParty/glslang/glslang/OSDependent/Windows/Debug/OSDependentd.lib)
	target_link_libraries(${APP_NAME} PRIVATE ${DILIGENT_ENGINE_BUILD_DIR}/DiligentCore/Primitives/Debug/Diligent-Primitives.lib)

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


set(SOURCE
    src/main.cpp
)

set(INCLUDE
    src/main.hpp
)

set(SHADERS
    assets/terrain.vsh
    assets/terrain.psh
    assets/terrain_wire.psh
    assets/terrain.gsh
    assets/terrain.hsh
    assets/terrain.dsh
    assets/structures.fxh
)

set(ASSETS
    assets/ps_height_1k.png
    assets/ps_texture_2k.png
)

add_sample_app("main" "./" "${SOURCE}" "${INCLUDE}" "${SHADERS}" "${ASSETS}")
# target_link_libraries(main PRIVATE Diligent-TextureLoader)

target_include_directories(main
PRIVATE
    src
    src/Terrain
	
	${DILIGENT_ENGINE_DIR}/DiligentCore/Common/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Platforms/Win32/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Platforms/Basic/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Platforms/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Primitives/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Graphics/GraphicsEngine/include
	${DILIGENT_ENGINE_DIR}/DiligentCore/Graphics/GraphicsEngine/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Graphics/GraphicsTools/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Graphics/GraphicsEngineD3D11/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Graphics/GraphicsEngineD3DBase/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Graphics/GraphicsEngineD3D12/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Graphics/GraphicsEngineOpenGL/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Graphics/GraphicsEngineVulkan/interface
	${DILIGENT_ENGINE_DIR}/DiligentCore/Graphics/GraphicsAccessories/interface
	
	${DILIGENT_ENGINE_DIR}/DiligentFX/.
	${DILIGENT_ENGINE_DIR}/DiligentFX/Components/interface
	${DILIGENT_ENGINE_DIR}/DiligentFX/PostProcess/EpipolarLightScattering/interface
	${DILIGENT_ENGINE_DIR}/DiligentFX/PostProcess/EpipolarLightScattering/../../Shaders/PostProcess/EpipolarLightScattering/public
	${DILIGENT_ENGINE_DIR}/DiligentFX/GLTF_PBR_Renderer/interface
	${DILIGENT_ENGINE_DIR}/DiligentFX/Components/interface
	${DILIGENT_ENGINE_DIR}/DiligentFX/PostProcess/EpipolarLightScattering/interface
    ${DILIGENT_ENGINE_DIR}/DiligentFX/Shaders/Common/public/
    ${DILIGENT_ENGINE_DIR}/DiligentFX/Shaders/PostProcess/EpipolarLightScattering/public
    ${DILIGENT_ENGINE_DIR}/DiligentFX/Shaders/PostProcess/ToneMapping/public
	
	${DILIGENT_ENGINE_DIR}/DiligentTools/AssetLoader/interface	
	${DILIGENT_ENGINE_DIR}/DiligentTools/Imgui/interface
	${DILIGENT_ENGINE_DIR}/DiligentTools/Imgui/../ThirdParty/imGuIZMO.quat
	${DILIGENT_ENGINE_DIR}/DiligentTools/Imgui/../ThirdParty/imgui	
	${DILIGENT_ENGINE_DIR}/DiligentTools/NativeApp/include
	${DILIGENT_ENGINE_DIR}/DiligentTools/NativeApp/include/Win32
	${DILIGENT_ENGINE_DIR}/DiligentTools/TextureLoader/interface
	${DILIGENT_ENGINE_DIR}/DiligentTools/ThirdParty/imgui
	${DILIGENT_ENGINE_DIR}/DiligentTools/ThirdParty/imGuIZMO.quat
	
	${DILIGENT_ENGINE_DIR}/DiligentSamples/SampleBase/include
	
)
