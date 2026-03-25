function build_website_documentation(options)
arguments
    options.rootDir = ".."
end

rootDir = char(java.io.File(char(options.rootDir)).getCanonicalPath());
buildFolder = fullfile(rootDir, "docs");
sourceFolder = fullfile(rootDir, "Documentation", "WebsiteDocumentation");

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
websiteFolder = "classes";
classes = {"Integrator", "IntegratorEulerMaruyama", "IntegratorWithDiffusivity", "IntegratorWithObstacles"};
classDocumentation = ClassDocumentation.empty(length(classes), 0);
for iName = 1:length(classes)
    classDocumentation(iName) = ClassDocumentation(classes{iName}, nav_order=iName, websiteRootURL=websiteRootURL, buildFolder=buildFolder, websiteFolder=websiteFolder, parent=classFolderName);
end
arrayfun(@(a) a.writeToFile(), classDocumentation)
trimTrailingWhitespaceInMarkdown(buildFolder)
end

function trimTrailingWhitespaceInMarkdown(rootFolder)
markdownFiles = dir(fullfile(rootFolder, "**", "*.md"));
for iFile = 1:numel(markdownFiles)
    filePath = fullfile(markdownFiles(iFile).folder, markdownFiles(iFile).name);
    fileText = fileread(filePath);
    trimmedText = regexprep(fileText, "[ \t]+(\r?\n)", "$1");
    if ~strcmp(fileText, trimmedText)
        fid = fopen(filePath, "w");
        assert(fid ~= -1, "Could not open markdown file for writing");
        fwrite(fid, trimmedText);
        fclose(fid);
    end
end
end
