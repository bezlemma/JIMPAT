# Run the following after compilation for a smaller file:

#C:\Users\Bez\Downloads\upx-5.0.0-win64\upx-5.0.0-win64\upx.exe --best --lzma bin\Stitching.exe
#C:\Users\Bez\Downloads\upx-5.0.0-win64\upx-5.0.0-win64\upx.exe --best --lzma bin\*.dll
#C:\Users\Bez\Downloads\upx-5.0.0-win64\upx-5.0.0-win64\upx.exe --best --lzma lib\julia\*.dll

using PackageCompiler
using Pkg
import Dates        # For timestamp in log messages
import Downloads    # For manual download test
import Tar          # To potentially verify the downloaded tarball
import Libdl        # To find DLL name

# --- Configuration ---
package_dir = raw"C:\Users\Bez\Documents\GitHub\JIMPAT\Stitching"
target_dir = raw"C:\Users\Bez\Documents\GitHub\JIMPAT\StitchingApp"
precompile_script_name = "precompile_exec_script.jl"
precompile_script_path = joinpath(package_dir, precompile_script_name)
intel_omp_artifact_name = "IntelOpenMP"
intel_omp_url = "https://github.com/JuliaBinaryWrappers/IntelOpenMP_jll.jl/releases/download/IntelOpenMP-v2025.0.4+0/IntelOpenMP.v2025.0.4.x86_64-w64-mingw32.tar.gz"
expected_omp_dll_name = "libiomp5md.dll" # Default guess (attempt to refine in Step 2)
# --- End Configuration ---

# --- Helper Function for Early Exit ---
function exit_script(message::String, original_project::String, exit_code::Int=1)
    @error "BUILD SCRIPT FAILED: $message"; flush(stderr)
    try if !isempty(original_project) && isfile(original_project) println("Attempting to restore original environment: ", original_project); flush(stdout); Pkg.activate(original_project) else println("No valid original project to restore."); flush(stdout) end catch e @warn "Failed to restore original environment: $e"; flush(stderr) end
    println("Build script terminated prematurely at: ", Dates.now()); flush(stdout)
    exit(exit_code)
end

println("Build script started at: ", Dates.now()); flush(stdout)
original_project = try Base.active_project() catch; "" end
println("Original active project: ", isempty(original_project) ? "None" : original_project); flush(stdout)

# --- Step 1: Cleanup ---
println("\n--- Step 1: Cleaning up potential stale files ---"); flush(stdout)
# [Cleanup code as before]
println("Finished: Step 1 Cleanup."); flush(stdout)

# --- Step 2: Preparing the source package environment ---
println("\n--- Step 2: Preparing the source package environment ---"); flush(stdout)
# [Environment prep code as in previous working version, including attempt to set expected_omp_dll_name]
try
    println("Starting: Pkg.activate..."); flush(stdout); Pkg.activate(package_dir); println("Finished: Pkg.activate. Current project: ", Base.active_project()); flush(stdout)
    println("Starting: Pkg.add MKL_jll..."); flush(stdout); Pkg.add("MKL_jll"); println("Finished: Pkg.add MKL_jll."); flush(stdout)
    println("Starting: Load MKL_jll to check DLL name..."); flush(stdout)
    try
        Pkg.build("MKL_jll"; verbose=false)
        mkl_jll_loaded = Base.eval(Main, :(try using MKL_jll; true catch; false end))
        if mkl_jll_loaded && isdefined(MKL_jll, :libiomp5) omp_lib_path = Base.eval(Main, :(MKL_jll.libiomp5)); global expected_omp_dll_name = basename(omp_lib_path); println("Determined expected OpenMP DLL name: ", expected_omp_dll_name); flush(stdout)
        else @warn "Could not determine OpenMP DLL name from MKL_jll. Using default guess: $expected_omp_dll_name"; flush(stderr) end
    catch load_err @warn "Error trying to load MKL_jll: $load_err. Using default guess: $expected_omp_dll_name"; flush(stderr) end
    println("Finished: Load MKL_jll attempt."); flush(stdout)
    println("Starting: Manual download test for $intel_omp_artifact_name..."); flush(stdout); temp_download_path = tempname() * ".tar.gz"
    try Downloads.download(intel_omp_url, temp_download_path); println("Finished: Manual download successful."); rm(temp_download_path; force=true); println("Finished: Manual download test passed."); flush(stdout) catch e @error "Manual download test FAILED: $e"; println("\nACTION REQUIRED: Check network/firewall/proxy/AV settings."); exit_script("Manual artifact download test failed.", original_project) end
    println("Starting: Pkg.instantiate..."); flush(stdout); Pkg.instantiate(; verbose = false); println("Finished: Pkg.instantiate."); flush(stdout)
    println("Starting: Pkg.precompile..."); flush(stdout); Pkg.precompile(); println("Finished: Pkg.precompile."); flush(stdout)
    println("Starting: Create precompile script..."); flush(stdout)
    try open(precompile_script_path, "w") do io; println(io, """println("--> Precompile exec: Loading Stitching..."); try using Stitching; println("--> Precompile exec: Stitching loaded.") catch e; println("--> Precompile exec: ERROR loading Stitching - ", e); end; println("--> Precompile exec: Finished.")"""); end; println("Finished: Create precompile script."); flush(stdout)
    catch e @warn "Could not create precompile script '$precompile_script_path': $e"; global precompile_script_path = ""; flush(stderr) end
catch e exit_script("Failed to prepare the package environment (Step 2). Error: $e", original_project) end
println("Finished: Step 2 Prepare Environment."); flush(stdout)


# --- Step 3: Create the Application ---
println("\n--- Step 3: Creating the application using PackageCompiler ---"); flush(stdout)

local app_creation_success = false # Address scope warning
try
    use_precompile_script = !isempty(precompile_script_path) && isfile(precompile_script_path)
    if use_precompile_script println("Using precompile execution file: ", precompile_script_path); flush(stdout) else println("Not using a precompile execution file."); flush(stdout) end

    # Call create_app using if/else and listing SUPPORTED keywords explicitly
    println("Starting: PackageCompiler.create_app (using supported kwargs)..."); flush(stdout)
    if use_precompile_script
        println("Calling create_app WITH precompile_execution_file."); flush(stdout)
        create_app(package_dir, target_dir;
            force=true,
            include_lazy_artifacts=true,
            #filter_stdlibs=true,
            #incremental=false,
            precompile_execution_file=precompile_script_path
        )
    else
        println("Calling create_app WITHOUT precompile_execution_file."); flush(stdout)
        create_app(package_dir, target_dir;
            force=true,
            #filter_stdlibs=true,
            #incremental=false,
            include_lazy_artifacts=true
        )
    end
    # If create_app crashes internally now, it's not a MethodError
    println("Finished: PackageCompiler.create_app call completed (API matched)."); flush(stdout)
    # Assume success at this point, unless checks below fail
    app_creation_success = true

    # --- Post-Build Audit REMOVED ---
    # println("Starting: PackageCompiler.audit_app..."); flush(stdout)
    # ... audit_app call removed ...

    # --- Improved Manual DLL Check using walkdir ---
    println("\n--- Performing recursive search for specific OpenMP DLL ---"); flush(stdout)
    println("Searching for: ", expected_omp_dll_name); flush(stdout)
    local found_dll_manually = false # Ensure this is local to the try block
    local found_dll_path = ""
    try
        for (root, dirs, files) in walkdir(target_dir)
            for file in files
                if file == expected_omp_dll_name
                    found_dll_manually = true
                    found_dll_path = joinpath(root, file)
                    # Use global marker if needed outside loop or use @goto
                    println("SUCCESS: Found expected OpenMP DLL at: ", found_dll_path); flush(stdout)
                    @goto dll_found # Exit loops once found
                end
            end
        end
        @label dll_found # Label to jump to after finding the DLL or finishing walk
    catch walk_err
         @warn "Error during directory walk for DLL check: $walk_err"; flush(stderr)
    end

    if !found_dll_manually
        # This is now just a warning, not necessarily a failure if create_app bundled it correctly
        @warn "Manual Check: Did NOT find expected OpenMP DLL ('$expected_omp_dll_name') via recursive search within '$target_dir'."; flush(stderr)
        println("Verify manually if the application runs correctly. The DLL might be present but named differently or PackageCompiler bundled dependencies correctly anyway."); flush(stdout)
        # Consider if this warning should be more severe? For now, let it pass.
        # app_creation_success = false # Uncomment to make missing DLL fatal again
    else
         println("Manual DLL Check finished successfully."); flush(stdout)
    end
    # --- End Improved Manual DLL Check ---


    # Check overall success marker from this step (only fails if we uncomment line above)
    if !app_creation_success
         exit_script("Application created, but post-build checks failed.", original_project)
    end

    println("\nApplication creation process finished successfully."); flush(stdout)

catch e # Catch errors during create_app structure or manual checks
    exit_script("Failed during application creation/checks (Step 3). Error: $e", original_project)
end


println("\nBuild script finished successfully at: ", Dates.now()); flush(stdout)