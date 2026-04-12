function build_website_documentation(options)
arguments
    options.rootDir = ".."
end

rootDir = char(java.io.File(char(options.rootDir)).getCanonicalPath());
buildFolder = fullfile(rootDir, "docs");
sourceFolder = fullfile(rootDir, "Documentation", "WebsiteDocumentation");

addPackageToPath(rootDir);

if isfolder(buildFolder)
    rmdir(buildFolder, "s");
end
copyfile(sourceFolder, buildFolder);

changelogPath = fullfile(rootDir, "CHANGELOG.md");
if isfile(changelogPath)
    header = "---" + newline + "layout: default" + newline + "title: Version History" + newline + "nav_order: 100" + newline + "---" + newline + newline;
    versionHistoryText = header + fileread(changelogPath);
    versionHistoryFilePath = fullfile(buildFolder, "version-history.md");
    fid = fopen(versionHistoryFilePath, "w");
    assert(fid ~= -1, "Could not open version-history.md for writing");
    fwrite(fid, versionHistoryText);
    fclose(fid);
end

evalin("base", "clear classes");
evalin("base", "rehash");

websiteRootURL = "advection-diffusion-models/";
classFolderName = "Class documentation";
classGroups = {
    struct( ...
        "parent", "Particle integrators", ...
        "grandparent", classFolderName, ...
        "websiteFolder", "classes/particle-integrators", ...
        "classes", {{"Integrator", "IntegratorEulerMaruyama", "IntegratorWithDiffusivity", "IntegratorWithObstacles", "AdvectionDiffusionIntegrator"}} ...
    )
    struct( ...
        "parent", "Kinematic models", ...
        "grandparent", classFolderName, ...
        "websiteFolder", "classes/kinematic-models", ...
        "classes", {{"KinematicModel", "StreamfunctionModel", "CylinderFlow", "DivergenceBox", "LinearVelocityField", "MeanderingJet", "NarrowEscapeProblem", "SimpleBox", "TranslatingGaussian"}} ...
    )
    struct( ...
        "parent", "Gridded streamfunction", ...
        "grandparent", "Estimators", ...
        "websiteFolder", "classes/estimators/gridded-streamfunction", ...
        "classes", {{"GriddedStreamfunction", "GriddedStreamfunctionBootstrap"}} ...
    )
    struct( ...
        "parent", "Density estimation", ...
        "grandparent", "Estimators", ...
        "websiteFolder", "classes/estimators/density-estimation", ...
        "classes", {{"KernelDensityEstimate"}} ...
    )
};

for iGroup = 1:numel(classGroups)
    group = classGroups{iGroup};
    classDocumentation = ClassDocumentation.empty(numel(group.classes), 0);
    for iName = 1:numel(group.classes)
        className = group.classes{iName};
        excludedMethodNames = string.empty(0, 1);
        switch className
            case {"KinematicModel", "StreamfunctionModel", "GriddedStreamfunction", "GriddedStreamfunctionBootstrap"}
                excludedMethodNames = string(className);
        end

        classDocumentation(iName) = ClassDocumentation( ...
            className, ...
            nav_order=iName, ...
            websiteRootURL=websiteRootURL, ...
            buildFolder=buildFolder, ...
            websiteFolder=group.websiteFolder, ...
            parent=group.parent, ...
            grandparent=group.grandparent, ...
            excludedMethodNames=excludedMethodNames);

        if ismember(className, ["GriddedStreamfunction", "GriddedStreamfunctionBootstrap"])
            classDocumentation(iName).addMethodDocumentation( ...
                annotatedWriteToFileDocumentation(className));
        end
    end
    arrayfun(@(a) a.writeToFile(), classDocumentation)
end

trimTrailingWhitespaceInMarkdown(fullfile(buildFolder, "classes", "estimators"))
end

function trimTrailingWhitespaceInMarkdown(rootFolder)
markdownFiles = dir(fullfile(rootFolder, "**", "*.md"));
for iFile = 1:numel(markdownFiles)
    filePath = fullfile(markdownFiles(iFile).folder, markdownFiles(iFile).name);
    fileText = fileread(filePath);
    trimmedText = regexprep(fileText, "[ \t]+(\r?\n)", "$1");
    trimmedText = regexprep(trimmedText, "(\r?\n){2,}$", newline);
    if ~strcmp(fileText, trimmedText)
        fid = fopen(filePath, "w");
        assert(fid ~= -1, "Could not open markdown file for writing");
        fwrite(fid, trimmedText);
        fclose(fid);
    end
end
end

function addPackageToPath(repoRoot)
metadataPath = fullfile(repoRoot, "resources", "mpackage.json");

if isfile(metadataPath)
    metadata = jsondecode(fileread(metadataPath));
    folderPaths = strings(length(metadata.folders), 1);
    for iFolder = 1:length(metadata.folders)
        folderPaths(iFolder) = string(fullfile(repoRoot, metadata.folders(iFolder).path));
    end
    addpath(folderPaths{:});
end
end

function metadata = annotatedWriteToFileDocumentation(className)
metadata = MethodDocumentation("writeToFile");
metadata.shortDescription = "Write this instance to a NetCDF restart file.";
metadata.topic = "Write to file";
metadata.declaration = "ncfile = writeToFile(self,path,properties=...,shouldOverwriteExisting=...,shouldAddRequiredProperties=...,attributes=...)";
metadata.parameters = [ ...
    struct("name", "path", "description", "path to the NetCDF restart file to write"), ...
    struct("name", "properties", "description", "optional property names to include in addition to the required restart state"), ...
    struct("name", "shouldOverwriteExisting", "description", "(optional) boolean indicating whether to overwrite an existing file at `path`. Default `false`."), ...
    struct("name", "shouldAddRequiredProperties", "description", "(optional) boolean indicating whether to include the required restart properties automatically. Default `true`."), ...
    struct("name", "attributes", "description", "(optional) dictionary of additional NetCDF attributes to attach to the root group")];
metadata.returns = struct("name", "ncfile", "description", "NetCDF file handle for the written restart file");
metadata.detailedDescription = char( ...
    "  Writes the canonical restart state for `" + className + ...
    "` to a NetCDF file so `fromFile` can reconstruct the saved object " + ...
    "later without rerunning the expensive fitting or bootstrap workflow." + ...
    newline + newline + ...
    "  Pass additional property names when you want to persist optional " + ...
    "state beyond the required restart payload.");
metadata.definingClassName = "CAAnnotatedClass";
metadata.addDeclaringClass(className);
metadata.functionType = FunctionType.instanceMethod;
end
