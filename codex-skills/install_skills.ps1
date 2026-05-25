$ErrorActionPreference = "Stop"

$skillNames = @(
    "binary-mask-stack",
    "nadph-phasor-trajectory",
    "tmrm-lgr5-nadph-analysis",
    "tmrm-nadph-analysis",
    "tmrm-nadph-percentile-plots"
)

$sourceRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$targetRoot = Join-Path $env:USERPROFILE ".codex\skills"

New-Item -ItemType Directory -Force -Path $targetRoot | Out-Null

foreach ($skillName in $skillNames) {
    $source = Join-Path $sourceRoot $skillName
    $target = Join-Path $targetRoot $skillName

    if (-not (Test-Path -LiteralPath $source -PathType Container)) {
        throw "Missing skill folder: $source"
    }

    Copy-Item -LiteralPath $source -Destination $target -Recurse -Force
    Write-Host "Installed $skillName"
}

Write-Host "Done. Restart Codex to reload the installed skills."
