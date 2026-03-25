function results = RunAllUnitTests()
testRoot = fileparts(mfilename('fullpath'));
originalFolder = pwd;
cleanup = onCleanup(@() cd(originalFolder));
cd(testRoot);

repoRoot = fileparts(testRoot);
workspaceRoot = fileparts(repoRoot);
dependencyRoot = fullfile(workspaceRoot,'spline-core');

addPackageToPath(repoRoot);

if ~isfolder(dependencyRoot)
    error('Expected dependency repo at %s.',dependencyRoot);
end
addPackageToPath(dependencyRoot);

import matlab.unittest.TestRunner
import matlab.unittest.TestSuite

runner = testrunner("textoutput");
suite = [ ...
    TestSuite.fromClass(?IntegratorUnitTests) ...
    TestSuite.fromClass(?IntegratorWithObstaclesUnitTests) ...
    TestSuite.fromClass(?LinearVelocityFieldUnitTests) ...
    TestSuite.fromClass(?LinearVelocityFieldEstimationUnitTests) ...
    ];

results = runner.run(suite);

if nargout == 0
    disp(table(results))
    if any(~[results.Passed])
        error('Some unit tests failed.');
    end
end
end

function addPackageToPath(repoRoot)
metadataPath = fullfile(repoRoot,'resources','mpackage.json');

if isfile(metadataPath)
    metadata = jsondecode(fileread(metadataPath));
    for iFolder = 1:length(metadata.folders)
        addpath(fullfile(repoRoot,metadata.folders(iFolder).path));
    end
end

addpath(repoRoot);
end
