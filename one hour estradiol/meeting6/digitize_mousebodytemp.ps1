param(
    [string]$InputPath = (Join-Path $PSScriptRoot 'references\mousebodytemp.png'),
    [string]$OutputPath = (Join-Path $PSScriptRoot 'references\mousebodytemp_digitized.csv')
)

Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Drawing

$panels = @(
    @{
        Name = 'nonpregnant_c'
        XStart = 116
        XStop = 662
        Slope = -0.025972565774679556
        Intercept = 46.484315268720486
    },
    @{
        Name = 'pregnant_control_c'
        XStart = 739
        XStop = 1285
        Slope = -0.01809942398514919
        Intercept = 43.65455488520652
    }
)

function Test-RedTracePixel {
    param([System.Drawing.Color]$Color)

    return (
        $Color.R -gt 110 -and
        ($Color.R - $Color.G) -gt 45 -and
        ($Color.R - $Color.B) -gt 45
    )
}

function Fill-Gaps {
    param([double[]]$Values)

    for ($index = 0; $index -lt $Values.Length; $index++) {
        if (-not [double]::IsNaN($Values[$index])) {
            continue
        }
        $left = $index - 1
        while ($left -ge 0 -and [double]::IsNaN($Values[$left])) { $left-- }
        $right = $index + 1
        while ($right -lt $Values.Length -and [double]::IsNaN($Values[$right])) { $right++ }

        if ($left -ge 0 -and $right -lt $Values.Length) {
            $fraction = ($index - $left) / ($right - $left)
            $Values[$index] = $Values[$left] + $fraction * ($Values[$right] - $Values[$left])
        } elseif ($left -ge 0) {
            $Values[$index] = $Values[$left]
        } elseif ($right -lt $Values.Length) {
            $Values[$index] = $Values[$right]
        }
    }
}

$resolvedInput = (Resolve-Path -LiteralPath $InputPath).Path
$bitmap = [System.Drawing.Bitmap]::FromFile($resolvedInput)
try {
    $traces = @{}
    foreach ($panel in $panels) {
        $slope = $panel.Slope
        $intercept = $panel.Intercept
        $values = [double[]]::new($panel.XStop - $panel.XStart + 1)

        for ($x = $panel.XStart; $x -le $panel.XStop; $x++) {
            $rowSum = 0.0
            $rowCount = 0
            for ($y = 250; $y -le 440; $y++) {
                if (Test-RedTracePixel -Color $bitmap.GetPixel($x, $y)) {
                    $rowSum += $y
                    $rowCount++
                }
            }
            $offset = $x - $panel.XStart
            $values[$offset] = if ($rowCount -gt 0) {
                $slope * ($rowSum / $rowCount) + $intercept
            } else {
                [double]::NaN
            }
        }
        Fill-Gaps -Values $values
        $traces[$panel.Name] = $values
    }
} finally {
    $bitmap.Dispose()
}

$count = $panels[0].XStop - $panels[0].XStart + 1
$records = for ($index = 0; $index -lt $count; $index++) {
    [pscustomobject]@{
        relative_day = 5.0 * $index / ($count - 1)
        nonpregnant_c = $traces.nonpregnant_c[$index]
        pregnant_control_c = $traces.pregnant_control_c[$index]
    }
}

$records | Export-Csv -LiteralPath $OutputPath -NoTypeInformation
[pscustomobject]@{
    Output = $OutputPath
    NonpregnantMin = ($records.nonpregnant_c | Measure-Object -Minimum).Minimum
    NonpregnantMax = ($records.nonpregnant_c | Measure-Object -Maximum).Maximum
    PregnantControlMin = ($records.pregnant_control_c | Measure-Object -Minimum).Minimum
    PregnantControlMax = ($records.pregnant_control_c | Measure-Object -Maximum).Maximum
}
